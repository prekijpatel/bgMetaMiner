import argparse
import json
import os
import random
import logging

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from torch.utils.data import Dataset, DataLoader
from transformers import AutoModel, AutoTokenizer, get_linear_schedule_with_warmup, logging as hf_logging

logging.basicConfig(filename="log_file.log", level=logging.DEBUG, format='%(asctime)s:%(levelname)s - %(message)s')
hf_logging.set_verbosity_error()


# The four levels the model predicts, and the three raw fields it reads.
label_columns = ["Manual_identified_host", "Manual_source_category", "Manual_source", "Manual_sample"]
input_columns = ["host", "host_disease", "isolation_source"]

# Friendlier names for the levels, used in tables and figure titles.
level_names = {"Manual_identified_host": "L1_host", "Manual_source_category": "L2_category",
               "Manual_source": "L3_source", "Manual_sample": "L4_sample"}

# Deep levels count double when we weigh how good a model is overall.
level_weights = {"Manual_identified_host": 1.0, "Manual_source_category": 1.0,
                 "Manual_source": 2.0, "Manual_sample": 2.0}

# The many ways "no value" shows up in raw metadata.
missing_tokens = {"", "na", "nan", "none", "missing", "not provided", "not applicable",
                  "not available", "not collected", "unknown", "n/a", "not determined", "not recorded"}

input_template = "[HOST] {host} [DIS] {disease} [SRC] {source}"
missing_placeholder = "[MISSING]"

seed = 42
default_model = "allenai/scibert_scivocab_uncased"
confidence_cutoffs = [0.9, 0.8, 0.7, 0.6, 0.5]


def set_seed(s:int=seed):
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(s)


def clean_field(value):
    # A single raw field -> lower-cased, or [MISSING] if it is empty/uninformative.
    if pd.isna(value):
        return missing_placeholder
    value = str(value).strip().lower()
    if value in missing_tokens:
        return missing_placeholder
    return value


def build_input_text(row:pd.Series):
    return input_template.format(host=clean_field(row.get("host")),
                                 disease=clean_field(row.get("host_disease")),
                                 source=clean_field(row.get("isolation_source")))


def load_and_prepare(tsv_path:str):
    logging.info(f"Reading the labelled metadata from {tsv_path}")
    df = pd.read_csv(tsv_path, sep="\t", low_memory=False)

    for col in label_columns + input_columns:
        if col not in df.columns:
            raise ValueError(f"The required column '{col}' is missing from {tsv_path}")

    # Every row needs a label at every level; blanks become "Unknown".
    for col in label_columns:
        df[col] = df[col].fillna("Unknown").astype(str).str.strip()
        df.loc[df[col] == "", col] = "Unknown"

    df["input_text"] = df.apply(build_input_text, axis=1)

    before = len(df)
    df = df.drop_duplicates(subset=["input_text"]).reset_index(drop=True)
    logging.info(f"Kept {len(df)} unique input strings out of {before} rows after removing exact duplicates.")
    return df


def fit_label_encoders(df:pd.DataFrame):
    # Fit on the full data so the class space is identical everywhere.
    encoders = {}
    for col in label_columns:
        enc = LabelEncoder()
        enc.fit(df[col].values)
        encoders[col] = enc
        rare = int((df[col].value_counts() < 5).sum())
        logging.info(f"{col}: {len(enc.classes_)} classes, {rare} of them with fewer than 5 examples.")
    return encoders


class MetaMinerDataset(Dataset):
    def __init__(self, df:pd.DataFrame, tokenizer, encoders:dict, max_length:int=64):
        self.texts = df["input_text"].tolist()
        self.tokenizer = tokenizer
        self.max_length = max_length
        self.labels = {col: encoders[col].transform(df[col].values) for col in label_columns}

    def __len__(self):
        return len(self.texts)

    def __getitem__(self, idx:int):
        encoded = self.tokenizer(self.texts[idx], max_length=self.max_length,
                                 padding="max_length", truncation=True, return_tensors="pt")
        item = {key: value.squeeze(0) for key, value in encoded.items()}
        for col in label_columns:
            item[col] = torch.tensor(self.labels[col][idx], dtype=torch.long)
        return item


class IsolationSourceClassifier(nn.Module):
    # One shared encoder feeding four independent linear heads, one per level.
    def __init__(self, model_name:str, classes_per_level:dict, dropout:float=0.1):
        super().__init__()
        self.encoder = AutoModel.from_pretrained(model_name, torch_dtype=torch.float32)
        hidden = self.encoder.config.hidden_size
        self.dropout = nn.Dropout(dropout)
        self.heads = nn.ModuleDict({col: nn.Linear(hidden, n) for col, n in classes_per_level.items()})

    def forward(self, input_ids, attention_mask, token_type_ids=None):
        kwargs = {"input_ids": input_ids, "attention_mask": attention_mask}
        if token_type_ids is not None:
            kwargs["token_type_ids"] = token_type_ids
        out = self.encoder(**kwargs)
        if getattr(out, "pooler_output", None) is not None:
            pooled = out.pooler_output
        else:
            pooled = out.last_hidden_state[:, 0, :]
        pooled = self.dropout(pooled)
        return {col: head(pooled) for col, head in self.heads.items()}


class FocalLoss(nn.Module):
    # Focuses the loss on hard/rare examples; an optional class-weight acts as alpha.
    def __init__(self, weight=None, gamma:float=2.0):
        super().__init__()
        self.weight = weight
        self.gamma = gamma

    def forward(self, logits, targets):
        log_prob = F.log_softmax(logits, dim=-1)
        true_log_prob = log_prob.gather(1, targets.unsqueeze(1)).squeeze(1)
        prob = true_log_prob.exp()
        loss = -((1.0 - prob) ** self.gamma) * true_log_prob
        if self.weight is not None:
            w = self.weight.gather(0, targets)
            return (loss * w).sum() / w.sum().clamp_min(1e-8)
        return loss.mean()


def class_weights(df:pd.DataFrame, col:str, n_classes:int, encoder:LabelEncoder, device):
    # Square-root inverse frequency, gentler than plain inverse, rescaled to mean 1.
    counts = np.zeros(n_classes, dtype=np.float64)
    for label, count in df[col].value_counts().items():
        counts[encoder.transform([label])[0]] = count
    weights = 1.0 / np.sqrt(counts + 1.0)
    weights = weights / weights.mean()
    return torch.tensor(weights, dtype=torch.float, device=device)


def train_model(train_df:pd.DataFrame, encoders:dict, model_name:str, args, device):
    logging.info(f"Training {model_name} on {len(train_df)} rows for {args.epochs} epochs (loss={args.loss}).")
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    classes_per_level = {col: len(encoders[col].classes_) for col in label_columns}

    dataset = MetaMinerDataset(train_df, tokenizer, encoders, args.max_length)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True,
                        num_workers=args.num_workers, pin_memory=True)

    model = IsolationSourceClassifier(model_name, classes_per_level).to(device)

    loss_fns = {}
    for col in label_columns:
        weight = None if args.no_class_weights else class_weights(train_df, col, classes_per_level[col], encoders[col], device)
        loss_fns[col] = FocalLoss(weight=weight, gamma=args.focal_gamma) if args.loss == "focal" else nn.CrossEntropyLoss(weight=weight)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=0.01)
    total_steps = args.epochs * max(1, len(loader))
    scheduler = get_linear_schedule_with_warmup(optimizer, int(0.1 * total_steps), total_steps)

    for epoch in range(1, args.epochs + 1):
        model.train()
        running_loss = 0.0
        for batch in loader:
            optimizer.zero_grad()
            input_ids = batch["input_ids"].to(device)
            attention_mask = batch["attention_mask"].to(device)
            token_type_ids = batch.get("token_type_ids")
            if token_type_ids is not None:
                token_type_ids = token_type_ids.to(device)

            logits = model(input_ids, attention_mask, token_type_ids)
            loss = sum(loss_fns[col](logits[col], batch[col].to(device)) for col in label_columns)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            scheduler.step()
            running_loss += loss.item()

        epoch_loss = running_loss / max(1, len(loader))
        print(f"    epoch {epoch}/{args.epochs}  loss={epoch_loss:.4f}")
        logging.debug(f"Epoch {epoch} finished with mean loss {epoch_loss:.4f}")

    return model, tokenizer


@torch.no_grad()
def predict_with_confidence(model, tokenizer, eval_df:pd.DataFrame, encoders:dict, args, device):
    # Returns, per level, the true labels, predicted labels and the top softmax confidence.
    logging.info(f"Scoring the held-out {len(eval_df)} rows.")
    model.eval()
    dataset = MetaMinerDataset(eval_df, tokenizer, encoders, args.max_length)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False,
                        num_workers=args.num_workers, pin_memory=True)

    truths = {col: [] for col in label_columns}
    preds = {col: [] for col in label_columns}
    confs = {col: [] for col in label_columns}

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        token_type_ids = batch.get("token_type_ids")
        if token_type_ids is not None:
            token_type_ids = token_type_ids.to(device)

        logits = model(input_ids, attention_mask, token_type_ids)
        for col in label_columns:
            probabilities = torch.softmax(logits[col], dim=-1)
            top_conf, top_idx = probabilities.max(dim=-1)
            preds[col].extend(top_idx.cpu().numpy().tolist())
            confs[col].extend(top_conf.cpu().numpy().tolist())
            truths[col].extend(batch[col].cpu().numpy().tolist())

    for col in label_columns:
        truths[col] = np.array(truths[col])
        preds[col] = np.array(preds[col])
        confs[col] = np.array(confs[col])
    return truths, preds, confs


def confidence_threshold_report(truths:dict, preds:dict, confs:dict, out_dir:str):
    # For every level and every cut-off, keep the predictions the model is at least
    # that confident about, and measure how good they are on that retained subset.
    rows = []
    total = len(truths[label_columns[0]])
    for col in label_columns:
        y_true, y_pred, c = truths[col], preds[col], confs[col]
        for cutoff in [0.0] + confidence_cutoffs:
            keep = c >= cutoff
            n_kept = int(keep.sum())
            if n_kept == 0:
                logging.warning(f"No predictions above confidence {cutoff} for {col}.")
                continue
            yt, yp = y_true[keep], y_pred[keep]
            rows.append({
                "level": level_names[col],
                "confidence_cutoff": cutoff,
                "coverage": round(n_kept / total, 4),
                "n_retained": n_kept,
                "accuracy": round(accuracy_score(yt, yp), 4),
                "recall_macro": round(recall_score(yt, yp, average="macro", zero_division=0), 4),
                "recall_weighted": round(recall_score(yt, yp, average="weighted", zero_division=0), 4),
                "f1_macro": round(f1_score(yt, yp, average="macro", zero_division=0), 4),
                "f1_weighted": round(f1_score(yt, yp, average="weighted", zero_division=0), 4),
                "precision_macro": round(precision_score(yt, yp, average="macro", zero_division=0), 4),
                "precision_weighted": round(precision_score(yt, yp, average="weighted", zero_division=0), 4),
            })

    report = pd.DataFrame(rows)
    path = os.path.join(out_dir, "holdout_confidence_metrics.tsv")
    report.to_csv(path, sep="\t", index=False)
    logging.info(f"Wrote the confidence-threshold report to {path}")

    print("\nHeld-out performance at each confidence cut-off (cut-off 0.00 = all predictions):")
    print(report.to_string(index=False))
    print(f"\nSaved: {path}")
    return report


def plot_class_distributions(df:pd.DataFrame, out_dir:str):
    # One panel per level: classes ranked by frequency, log y-axis to show the long tail.
    logging.info("Drawing the per-level class-frequency distributions.")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, col in zip(axes.ravel(), label_columns):
        counts = df[col].value_counts().sort_values(ascending=False)
        ax.bar(range(1, len(counts) + 1), counts.values, width=1.0, color="#4C72B0")
        ax.set_yscale("log")
        ax.set_title(f"{level_names[col]}  ({len(counts)} classes)")
        ax.set_xlabel("class rank (most to least frequent)")
        ax.set_ylabel("number of records (log)")
        rare = int((counts < 5).sum())
        ax.text(0.97, 0.95, f"{rare} classes < 5 records", transform=ax.transAxes,
                ha="right", va="top", fontsize=9, color="grey")
        # The underlying counts, in case the numbers are needed for the manuscript.
        counts.rename_axis("class").reset_index(name="count").to_csv(
            os.path.join(out_dir, f"class_distribution_{level_names[col]}.tsv"), sep="\t", index=False)

    fig.suptitle("Class-frequency distribution per classification level", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out_dir, f"class_distribution.{ext}"), dpi=300)
    plt.close(fig)
    logging.info(f"Saved the class-distribution figure and per-level count tables to {out_dir}")
    print(f"Saved: {os.path.join(out_dir, 'class_distribution.png')} (and .pdf, plus per-level count TSVs)")


def save_deployable_model(model, tokenizer, encoders:dict, model_name:str, args, save_dir:str):
    # Everything the in-tool inference wrapper needs to run offline.
    os.makedirs(save_dir, exist_ok=True)
    torch.save(model.state_dict(), os.path.join(save_dir, "pytorch_model.bin"))
    tokenizer.save_pretrained(os.path.join(save_dir, "tokenizer"))
    model.encoder.config.save_pretrained(os.path.join(save_dir, "encoder"))

    with open(os.path.join(save_dir, "label_encoders.json"), "w", encoding="utf-8") as f:
        json.dump({col: encoders[col].classes_.tolist() for col in label_columns}, f, indent=2, ensure_ascii=False)

    config = {"base_model_name": model_name, "label_columns": label_columns, "input_columns": input_columns,
              "n_classes_per_head": {col: len(encoders[col].classes_) for col in label_columns},
              "max_length": args.max_length, "input_template": input_template,
              "missing_placeholder": missing_placeholder, "missing_tokens": sorted(missing_tokens)}
    with open(os.path.join(save_dir, "config.json"), "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)
    logging.info(f"Saved the deployable model to {save_dir}")
    print(f"Saved deployable model: {save_dir}")


def main():
    parser = argparse.ArgumentParser(description="Hold-out evaluation of the isolation-source NLP classifier")
    parser.add_argument("--data", required=True, help="Path to the labelled TSV")
    parser.add_argument("--output", default="holdout_results", help="Where to write reports and figures")
    parser.add_argument("--model_name", default=default_model)
    parser.add_argument("--test_size", type=float, default=0.10, help="Held-out fraction (0.10 = 90/10 split)")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--lr", type=float, default=2e-5)
    parser.add_argument("--max_length", type=int, default=64)
    parser.add_argument("--num_workers", type=int, default=6)
    parser.add_argument("--loss", choices=["ce", "focal"], default="focal")
    parser.add_argument("--focal_gamma", type=float, default=2.0)
    parser.add_argument("--no_class_weights", action="store_true")
    parser.add_argument("--save_final_to", default=None,
                        help="If set, also retrain on the FULL dataset and save the deployable model here")
    args = parser.parse_args()

    set_seed()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}  |  loss: {args.loss}")
    logging.info(f"Starting hold-out evaluation on {device} with loss={args.loss}.")
    if device.type == "cpu":
        print("NOTE: no GPU found, this will be slow.")

    os.makedirs(args.output, exist_ok=True)

    df = load_and_prepare(args.data)
    encoders = fit_label_encoders(df)

    # The long-tail figure is about the taxonomy itself, so it uses the full data.
    plot_class_distributions(df, args.output)

    # 90/10 split, stratified on the host level so every host is represented on both sides.
    train_df, test_df = train_test_split(df, test_size=args.test_size, random_state=seed,
                                         stratify=df["Manual_identified_host"])
    train_df = train_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)
    print(f"Split: train={len(train_df)}  held-out test={len(test_df)}")
    logging.info(f"90/10 split gave {len(train_df)} training and {len(test_df)} held-out rows.")

    model, tokenizer = train_model(train_df, encoders, args.model_name, args, device)
    truths, preds, confs = predict_with_confidence(model, tokenizer, test_df, encoders, args, device)
    confidence_threshold_report(truths, preds, confs, args.output)

    # Optionally produce the model that actually ships, trained on everything.
    if args.save_final_to:
        print("\nRetraining on the full dataset for the deployable model...")
        logging.info("Retraining on the full dataset for deployment.")
        full_model, full_tokenizer = train_model(df, encoders, args.model_name, args, device)
        save_deployable_model(full_model, full_tokenizer, encoders, args.model_name, args, args.save_final_to)

    print(f"\nDone. Everything is in {os.path.abspath(args.output)}")
    logging.info("Hold-out evaluation finished.")


if __name__ == "__main__":
    main()
