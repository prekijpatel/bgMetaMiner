
import argparse
import contextlib
import json
import os
import random
import warnings
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import accuracy_score, f1_score, classification_report
from sklearn.model_selection import StratifiedKFold, train_test_split
from sklearn.preprocessing import LabelEncoder
from torch.utils.data import DataLoader, Dataset
from transformers import AutoModel, AutoTokenizer, get_linear_schedule_with_warmup, logging as hf_logging

logging.basicConfig(filename="log_file.log", level=logging.DEBUG, format='%(asctime)s:%(levelname)s - %(message)s')
warnings.filterwarnings("ignore")
hf_logging.set_verbosity_error()
# Keeps the HF tokenizer quiet when DataLoader workers fork.
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

# MLflow is optional, so the script still runs if it is not installed or is
# switched off with --no_mlflow.
try:
    import mlflow
    mlflow_available = True
except ImportError:
    mlflow_available = False

mlflow_enabled = False  # decided in main() from the args and availability


# The candidate encoders, the four target levels, and the three raw input fields.
models_to_benchmark = {
    "deberta-v3-small": "microsoft/deberta-v3-small",
    "pubmedbert": "microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext",
    "scibert": "allenai/scibert_scivocab_uncased",
}

label_columns = ["Manual_identified_host", "Manual_source_category", "Manual_source_y", "Manual_sample"]
input_columns = ["host", "host_disease", "isolation_source"]

# Friendlier level names used in the tables, figures and MLflow metric keys.
level_names = {"Manual_identified_host": "L1_host", "Manual_source_category": "L2_category",
               "Manual_source_y": "L3_source", "Manual_sample": "L4_sample"}

# The composite selection score is a weighted mean of the per-level macro-F1s.
# Deep levels count double, since L3/L4 are the hard, useful part.
level_weights = {"Manual_identified_host": 1.0, "Manual_source_category": 1.0,
                 "Manual_source_y": 2.0, "Manual_sample": 2.0}

# The many spellings of "no value" in raw metadata.
missing_tokens = {"", "na", "nan", "none", "missing", "not provided", "not applicable",
                  "not available", "not collected", "unknown", "n/a", "not determined", "not recorded"}

input_template = "[HOST] {host} [DIS] {disease} [SRC] {source}"
missing_placeholder = "[MISSING]"
seed = 42

# Filled in main() so class-weighting can map label strings back to their indices.
encoder_lookup = {}


@contextlib.contextmanager
def mlflow_run(run_name=None, nested:bool=False):
    if mlflow_enabled:
        # A left-over active run (e.g. after a Ctrl-C) would block a new top-level
        # run, so close it first.
        if not nested and mlflow.active_run() is not None:
            mlflow.end_run()
        with mlflow.start_run(run_name=run_name, nested=nested) as run:
            yield run
    else:
        yield None


def mlflow_log_params(params:dict):
    if mlflow_enabled:
        mlflow.log_params(params)


def mlflow_log_metrics(metrics:dict, step=None):
    if mlflow_enabled:
        mlflow.log_metrics(metrics, step=step)


def mlflow_log_artifact(path):
    if mlflow_enabled and Path(path).exists():
        mlflow.log_artifact(str(path))


def set_seed(s:int=seed):
    random.seed(s)
    np.random.seed(s)
    torch.manual_seed(s)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(s)


def clean_field(value):
    # A single raw field -> lower-cased, or [MISSING] if empty/uninformative.
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
    print(f"Loaded {before} rows -> {len(df)} unique input strings")
    logging.info(f"Kept {len(df)} unique input strings out of {before} rows after removing exact duplicates.")
    return df


def fit_label_encoders(df:pd.DataFrame):
    # Fit on the full data so the class space is identical across folds, the
    # held-out test set and the final model.
    encoders = {}
    print("\nLabel distribution per level:")
    for col in label_columns:
        enc = LabelEncoder()
        enc.fit(df[col].values)
        encoders[col] = enc
        rare = int((df[col].value_counts() < 5).sum())
        print(f"  {col:30s} {len(enc.classes_):3d} classes, {rare} have <5 examples")
        logging.info(f"{col}: {len(enc.classes_)} classes, {rare} with fewer than 5 examples.")
    return encoders


def stratified_split(df:pd.DataFrame, stratify_col:str, test_size:float=0.15, random_state:int=seed):
    # Falls back to a plain random split if some class is too small to stratify.
    if df[stratify_col].value_counts().min() < 2:
        logging.warning(f"Some classes in {stratify_col} have <2 examples; using a random split instead.")
        return train_test_split(df, test_size=test_size, random_state=random_state)
    return train_test_split(df, test_size=test_size, stratify=df[stratify_col], random_state=random_state)


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
        # torch_dtype=float32 pins the encoder to 32-bit; newer transformers may
        # otherwise load fp16 weights that clash with the fp32 linear heads.
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
    # Focuses the loss on hard/rare examples; an optional class-weight acts as
    # alpha and composes with the sqrt-inverse-frequency weights. Weighted
    # reduction matches nn.CrossEntropyLoss(weight=...): sum(w*loss)/sum(w).
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


def compute_class_weights(df:pd.DataFrame, col:str, n_classes:int, device):
    # Square-root inverse frequency, aligned to the LabelEncoder order and
    # rescaled to a mean of 1.
    counts = np.zeros(n_classes, dtype=np.float64)
    for label, count in df[col].value_counts().items():
        counts[encoder_lookup[col].transform([label])[0]] = count
    weights = 1.0 / np.sqrt(counts + 1.0)
    weights = weights / weights.mean()
    return torch.tensor(weights, dtype=torch.float, device=device)


def weighted_macro_f1(per_level:dict):
    # Deep-level-weighted mean of the per-level macro-F1s; drives selection and
    # early stopping.
    numerator = sum(level_weights[c] * per_level[c]["macro_f1"] for c in label_columns)
    return numerator / sum(level_weights[c] for c in label_columns)


def make_loaders(train_df:pd.DataFrame, eval_df:pd.DataFrame, tokenizer, encoders:dict, args):
    train_ds = MetaMinerDataset(train_df, tokenizer, encoders, args.max_length)
    eval_ds = MetaMinerDataset(eval_df, tokenizer, encoders, args.max_length)
    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True,
                              num_workers=args.num_workers, pin_memory=True)
    eval_loader = DataLoader(eval_ds, batch_size=args.batch_size, shuffle=False,
                             num_workers=args.num_workers, pin_memory=True)
    return train_loader, eval_loader


def make_loss_fns(train_df:pd.DataFrame, n_classes:dict, device, args):
    # One loss per head: cross-entropy or focal, each with optional class weights.
    use_weights = not args.no_class_weights
    loss_fns = {}
    for col in label_columns:
        weight = compute_class_weights(train_df, col, n_classes[col], device) if use_weights else None
        if args.loss == "focal":
            loss_fns[col] = FocalLoss(weight=weight, gamma=args.focal_gamma)
        else:
            loss_fns[col] = nn.CrossEntropyLoss(weight=weight)
    return loss_fns


def make_optimizer_scheduler(model, n_steps_per_epoch:int, epochs:int, lr:float):
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=0.01)
    total_steps = epochs * max(1, n_steps_per_epoch)
    scheduler = get_linear_schedule_with_warmup(optimizer, int(0.1 * total_steps), total_steps)
    return optimizer, scheduler


def train_one_epoch(model, loader, optimizer, scheduler, loss_fns:dict, device):
    model.train()
    total_loss = 0.0
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
        total_loss += loss.item()
    return total_loss / max(1, len(loader))


@torch.no_grad()
def evaluate(model, loader, encoders:dict, device):
    model.eval()
    preds = {col: [] for col in label_columns}
    targets = {col: [] for col in label_columns}

    for batch in loader:
        input_ids = batch["input_ids"].to(device)
        attention_mask = batch["attention_mask"].to(device)
        token_type_ids = batch.get("token_type_ids")
        if token_type_ids is not None:
            token_type_ids = token_type_ids.to(device)

        logits = model(input_ids, attention_mask, token_type_ids)
        for col in label_columns:
            preds[col].extend(logits[col].argmax(dim=-1).cpu().numpy().tolist())
            targets[col].extend(batch[col].cpu().numpy().tolist())

    per_level = {}
    for col in label_columns:
        y_true, y_pred = np.array(targets[col]), np.array(preds[col])
        per_level[col] = {
            "accuracy": float(accuracy_score(y_true, y_pred)),
            "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
            "weighted_f1": float(f1_score(y_true, y_pred, average="weighted", zero_division=0)),
            "n_classes": int(len(encoders[col].classes_)),
        }

    # Full-path accuracy: all four levels correct at once.
    n = len(targets[label_columns[0]])
    all_correct = sum(all(preds[col][i] == targets[col][i] for col in label_columns) for i in range(n))
    full_path = all_correct / max(1, n)
    return per_level, full_path, preds, targets


def per_class_report(preds:dict, targets:dict, encoders:dict):
    # Per-class precision/recall/F1/support for each level, as one DataFrame.
    rows = []
    for col in label_columns:
        y_true, y_pred = np.array(targets[col]), np.array(preds[col])
        labels_present = sorted(set(y_true.tolist()) | set(y_pred.tolist()))
        report = classification_report(y_true, y_pred, labels=labels_present,
                                       target_names=[encoders[col].classes_[i] for i in labels_present],
                                       output_dict=True, zero_division=0)
        for cls_name, metrics in report.items():
            if not isinstance(metrics, dict):
                continue
            rows.append({"level": col, "class": cls_name,
                         "precision": metrics.get("precision", 0.0),
                         "recall": metrics.get("recall", 0.0),
                         "f1": metrics.get("f1-score", 0.0),
                         "support": metrics.get("support", 0)})
    return pd.DataFrame(rows)


def log_epoch_metrics(per_level:dict, full_path:float, score:float, train_loss:float, epoch:int):
    metrics = {"train_loss": train_loss, "val_full_path_acc": full_path, "val_selection_score": score}
    for col in label_columns:
        short = level_names[col]
        metrics[f"val_{short}_acc"] = per_level[col]["accuracy"]
        metrics[f"val_{short}_macroF1"] = per_level[col]["macro_f1"]
        metrics[f"val_{short}_weightedF1"] = per_level[col]["weighted_f1"]
    mlflow_log_metrics(metrics, step=epoch)


def train_with_early_stopping(model, train_loader, val_loader, loss_fns:dict, encoders:dict,
                              args, device, restore_best:bool=True):
    # Trains up to args.max_epochs, stopping once the selection score has not
    # improved for args.patience epochs. Keeps the best epoch's numbers (and,
    # if asked, its weights).
    optimizer, scheduler = make_optimizer_scheduler(model, len(train_loader), args.max_epochs, args.lr)

    best = {"selection_score": -1.0, "best_epoch": 0, "full_path": 0.0,
            "per_level": None, "preds": None, "targets": None}
    best_state = None
    epochs_no_improve = 0

    for epoch in range(1, args.max_epochs + 1):
        train_loss = train_one_epoch(model, train_loader, optimizer, scheduler, loss_fns, device)
        per_level, path_acc, preds, targets = evaluate(model, val_loader, encoders, device)
        score = weighted_macro_f1(per_level)
        log_epoch_metrics(per_level, path_acc, score, train_loss, epoch)

        improved = score > best["selection_score"] + args.min_delta
        marker = "  *best" if improved else ""
        l3, l4 = per_level["Manual_source_y"], per_level["Manual_sample"]  # shown as macro/weighted
        print(f"    epoch {epoch:2d}/{args.max_epochs}  loss={train_loss:.4f}  "
              f"score={score:.4f}  full_path={path_acc:.4f}"
              f"  L3[m/w]={l3['macro_f1']:.3f}/{l3['weighted_f1']:.3f}"
              f"  L4[m/w]={l4['macro_f1']:.3f}/{l4['weighted_f1']:.3f}{marker}")

        if improved:
            best.update({"selection_score": score, "best_epoch": epoch, "full_path": path_acc,
                         "per_level": per_level, "preds": preds, "targets": targets})
            if restore_best:
                best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= args.patience:
                print(f"    early stop at epoch {epoch} (no improvement for {args.patience} "
                      f"epochs; best was epoch {best['best_epoch']})")
                logging.info(f"Early stopped at epoch {epoch}; best epoch was {best['best_epoch']}.")
                break

    if restore_best and best_state is not None:
        model.load_state_dict(best_state)
    return best


def aggregate_folds(fold_bests:list):
    # Mean and std of every metric across the folds.
    def mean_std(values):
        return float(np.mean(values)), float(np.std(values))

    agg = {"n_folds": len(fold_bests)}
    agg["best_epoch_mean"], agg["best_epoch_std"] = mean_std([b["best_epoch"] for b in fold_bests])
    agg["selection_score_mean"], agg["selection_score_std"] = mean_std([b["selection_score"] for b in fold_bests])
    agg["full_path_mean"], agg["full_path_std"] = mean_std([b["full_path"] for b in fold_bests])

    agg["levels"] = {}
    for col in label_columns:
        acc_m, acc_s = mean_std([b["per_level"][col]["accuracy"] for b in fold_bests])
        f1_m, f1_s = mean_std([b["per_level"][col]["macro_f1"] for b in fold_bests])
        w_m, w_s = mean_std([b["per_level"][col]["weighted_f1"] for b in fold_bests])
        agg["levels"][col] = {"accuracy_mean": acc_m, "accuracy_std": acc_s,
                              "macro_f1_mean": f1_m, "macro_f1_std": f1_s,
                              "weighted_f1_mean": w_m, "weighted_f1_std": w_s,
                              "n_classes": fold_bests[0]["per_level"][col]["n_classes"]}
    return agg


def log_fold_aggregate(agg:dict):
    metrics = {"cv_selection_score_mean": agg["selection_score_mean"],
               "cv_selection_score_std": agg["selection_score_std"],
               "cv_full_path_mean": agg["full_path_mean"],
               "cv_full_path_std": agg["full_path_std"],
               "cv_best_epoch_mean": agg["best_epoch_mean"]}
    for col in label_columns:
        short = level_names[col]
        lv = agg["levels"][col]
        metrics[f"cv_{short}_macroF1_mean"] = lv["macro_f1_mean"]
        metrics[f"cv_{short}_macroF1_std"] = lv["macro_f1_std"]
        metrics[f"cv_{short}_weightedF1_mean"] = lv["weighted_f1_mean"]
        metrics[f"cv_{short}_weightedF1_std"] = lv["weighted_f1_std"]
        metrics[f"cv_{short}_acc_mean"] = lv["accuracy_mean"]
    mlflow_log_metrics(metrics)


def print_cv_summary(model_key:str, agg:dict):
    print(f"\n  --- {model_key} cross-validation summary (mean +/- std over {agg['n_folds']} folds) ---")
    print(f"    selection score : {agg['selection_score_mean']:.4f} +/- {agg['selection_score_std']:.4f}")
    print(f"    full-path acc   : {agg['full_path_mean']:.4f} +/- {agg['full_path_std']:.4f}")
    print(f"    best epoch      : {agg['best_epoch_mean']:.1f} +/- {agg['best_epoch_std']:.1f}")
    for col in label_columns:
        lv = agg["levels"][col]
        print(f"    {level_names[col]:12s} acc={lv['accuracy_mean']:.4f}  "
              f"macroF1={lv['macro_f1_mean']:.4f}+/-{lv['macro_f1_std']:.4f}  "
              f"weightedF1={lv['weighted_f1_mean']:.4f}+/-{lv['weighted_f1_std']:.4f}  "
              f"({lv['n_classes']} classes)")


def run_cross_validation(model_key:str, model_name:str, dev_df:pd.DataFrame, encoders:dict, args, device):
    print(f"\n{'=' * 78}\n  CROSS-VALIDATION: {model_key}  ({model_name})  --  {args.k_folds} folds\n{'=' * 78}")
    logging.info(f"Cross-validating {model_key} ({model_name}) over {args.k_folds} folds.")

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    n_classes = {col: len(encoders[col].classes_) for col in label_columns}

    skf = StratifiedKFold(n_splits=args.k_folds, shuffle=True, random_state=seed)
    strat = dev_df["Manual_identified_host"].values  # stratify on the host level

    fold_bests = []
    with mlflow_run(run_name=f"cv_{model_key}_{args.loss}"):
        mlflow_log_params({"model_key": model_key, "model_name": model_name, "mode": "cross_validation",
                           "k_folds": args.k_folds, "max_epochs": args.max_epochs, "patience": args.patience,
                           "min_delta": args.min_delta, "lr": args.lr, "batch_size": args.batch_size,
                           "max_length": args.max_length, "loss": args.loss, "focal_gamma": args.focal_gamma,
                           "class_weights": not args.no_class_weights, "level_weights": json.dumps(level_weights),
                           "n_dev": len(dev_df),
                           **{f"n_classes_{level_names[c]}": n_classes[c] for c in label_columns}})

        for fold, (train_idx, val_idx) in enumerate(skf.split(dev_df, strat), 1):
            train_df = dev_df.iloc[train_idx].reset_index(drop=True)
            val_df = dev_df.iloc[val_idx].reset_index(drop=True)
            print(f"\n  Fold {fold}/{args.k_folds}  train={len(train_df)}  val={len(val_df)}")

            train_loader, val_loader = make_loaders(train_df, val_df, tokenizer, encoders, args)
            model = IsolationSourceClassifier(model_name, n_classes).to(device)
            loss_fns = make_loss_fns(train_df, n_classes, device, args)

            with mlflow_run(run_name=f"{model_key}_{args.loss}_fold{fold}", nested=True):
                mlflow_log_params({"fold": fold, "n_train": len(train_df), "n_val": len(val_df)})
                best = train_with_early_stopping(model, train_loader, val_loader, loss_fns,
                                                 encoders, args, device, restore_best=False)
                fold_summary = {"fold_best_epoch": best["best_epoch"],
                                "fold_selection_score": best["selection_score"],
                                "fold_full_path": best["full_path"]}
                for col in label_columns:
                    fold_summary[f"fold_{level_names[col]}_macroF1"] = best["per_level"][col]["macro_f1"]
                    fold_summary[f"fold_{level_names[col]}_weightedF1"] = best["per_level"][col]["weighted_f1"]
                mlflow_log_metrics(fold_summary)

            fold_bests.append(best)
            del model  # free the GPU between folds
            if device.type == "cuda":
                torch.cuda.empty_cache()

        agg = aggregate_folds(fold_bests)
        log_fold_aggregate(agg)

    print_cv_summary(model_key, agg)
    return agg


def retrain_and_test(model_key:str, dev_df:pd.DataFrame, test_df:pd.DataFrame, encoders:dict,
                     args, device, epochs:int, out_dir:Path):
    # Retrain the winner on the full dev set for a fixed epoch count, then score
    # it once on the untouched test set. These are the honest, reportable numbers.
    model_name = models_to_benchmark[model_key]
    print(f"\n{'=' * 78}\n  FINAL TEST EVALUATION: {model_key} "
          f"(retrain on dev for {epochs} epochs, eval on held-out test)\n{'=' * 78}")
    logging.info(f"Final test evaluation for {model_key} on the held-out set.")

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    n_classes = {col: len(encoders[col].classes_) for col in label_columns}
    train_loader, test_loader = make_loaders(dev_df, test_df, tokenizer, encoders, args)

    model = IsolationSourceClassifier(model_name, n_classes).to(device)
    loss_fns = make_loss_fns(dev_df, n_classes, device, args)
    optimizer, scheduler = make_optimizer_scheduler(model, len(train_loader), epochs, args.lr)

    for epoch in range(1, epochs + 1):
        train_loss = train_one_epoch(model, train_loader, optimizer, scheduler, loss_fns, device)
        print(f"    epoch {epoch}/{epochs}  loss={train_loss:.4f}")
        mlflow_log_metrics({"test_retrain_loss": train_loss}, step=epoch)

    per_level, path_acc, preds, targets = evaluate(model, test_loader, encoders, device)

    print(f"\n  HELD-OUT TEST RESULTS ({model_key}):")
    print(f"    full-path accuracy : {path_acc:.4f}")
    test_metrics = {"test_full_path_acc": path_acc, "test_selection_score": weighted_macro_f1(per_level)}
    for col in label_columns:
        m = per_level[col]
        short = level_names[col]
        print(f"    {short:12s} acc={m['accuracy']:.4f}  macroF1={m['macro_f1']:.4f}  "
              f"weightedF1={m['weighted_f1']:.4f}  ({m['n_classes']} classes)")
        test_metrics[f"test_{short}_acc"] = m["accuracy"]
        test_metrics[f"test_{short}_macroF1"] = m["macro_f1"]
        test_metrics[f"test_{short}_weightedF1"] = m["weighted_f1"]
    mlflow_log_metrics(test_metrics)

    pc_df = per_class_report(preds, targets, encoders)
    pc_path = out_dir / f"per_class_test_{model_key}.tsv"
    pc_df.to_csv(pc_path, sep="\t", index=False)
    mlflow_log_artifact(pc_path)
    print(f"    per-class test report -> {pc_path}")
    return per_level, path_acc


def build_comparison(all_aggs:dict, out_dir:Path):
    rows = []
    for key, agg in all_aggs.items():
        row = {"model": key, "best_epoch_mean": round(agg["best_epoch_mean"], 1),
               "selection_score": f"{agg['selection_score_mean']:.4f}+/-{agg['selection_score_std']:.4f}",
               "full_path": f"{agg['full_path_mean']:.4f}+/-{agg['full_path_std']:.4f}"}
        for col in label_columns:
            short = level_names[col]
            lv = agg["levels"][col]
            row[f"{short}_acc"] = f"{lv['accuracy_mean']:.4f}"
            row[f"{short}_macroF1"] = f"{lv['macro_f1_mean']:.4f}+/-{lv['macro_f1_std']:.4f}"
            row[f"{short}_weightedF1"] = f"{lv['weighted_f1_mean']:.4f}+/-{lv['weighted_f1_std']:.4f}"
        rows.append(row)

    comparison = pd.DataFrame(rows)
    print("\n" + "=" * 78)
    print("  CROSS-VALIDATION COMPARISON (mean +/- std over folds, all levels)")
    print("=" * 78)
    print(comparison.to_string(index=False))

    path = out_dir / "cv_comparison.tsv"
    comparison.to_csv(path, sep="\t", index=False)
    mlflow_log_artifact(path)
    return comparison


def save_artifacts(save_dir, model, tokenizer, encoders:dict, model_name:str, args):
    # Everything the in-tool inference wrapper needs to run offline.
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    torch.save(model.state_dict(), save_dir / "pytorch_model.bin")
    tokenizer.save_pretrained(save_dir / "tokenizer")
    model.encoder.config.save_pretrained(save_dir / "encoder")

    with open(save_dir / "label_encoders.json", "w", encoding="utf-8") as f:
        json.dump({col: encoders[col].classes_.tolist() for col in label_columns}, f, indent=2, ensure_ascii=False)

    config = {"base_model_name": model_name, "label_columns": label_columns, "input_columns": input_columns,
              "n_classes_per_head": {col: len(encoders[col].classes_) for col in label_columns},
              "max_length": args.max_length, "input_template": input_template,
              "missing_placeholder": missing_placeholder, "missing_tokens": sorted(missing_tokens)}
    with open(save_dir / "config.json", "w", encoding="utf-8") as f:
        json.dump(config, f, indent=2)

    print(f"  [SAVED] deployable artifacts -> {save_dir.resolve()}")
    logging.info(f"Saved the deployable model to {save_dir.resolve()}")
    mlflow_log_artifact(save_dir / "config.json")
    mlflow_log_artifact(save_dir / "label_encoders.json")


def train_final_model(model_key:str, df:pd.DataFrame, encoders:dict, args, device, save_dir):
    # Train the chosen encoder on the FULL dataset for a fixed epoch count, save.
    model_name = models_to_benchmark[model_key]
    epochs = args.final_epochs if args.final_epochs else args.max_epochs
    print(f"\n{'=' * 78}\n  FINAL TRAINING: {model_key} ({model_name}) "
          f"on {len(df)} rows for {epochs} epochs\n{'=' * 78}")
    logging.info(f"Training the final {model_key} model on {len(df)} rows for {epochs} epochs.")

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    n_classes = {col: len(encoders[col].classes_) for col in label_columns}

    dataset = MetaMinerDataset(df, tokenizer, encoders, args.max_length)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True,
                        num_workers=args.num_workers, pin_memory=True)

    model = IsolationSourceClassifier(model_name, n_classes).to(device)
    loss_fns = make_loss_fns(df, n_classes, device, args)
    optimizer, scheduler = make_optimizer_scheduler(model, len(loader), epochs, args.lr)

    mlflow_log_params({"model_key": model_key, "model_name": model_name, "mode": "train_final",
                       "final_epochs": epochs, "n_rows": len(df), "lr": args.lr,
                       "batch_size": args.batch_size, "max_length": args.max_length,
                       "loss": args.loss, "focal_gamma": args.focal_gamma,
                       "class_weights": not args.no_class_weights})

    for epoch in range(1, epochs + 1):
        train_loss = train_one_epoch(model, loader, optimizer, scheduler, loss_fns, device)
        print(f"    epoch {epoch}/{epochs}  loss={train_loss:.4f}")
        mlflow_log_metrics({"final_train_loss": train_loss}, step=epoch)

    save_artifacts(save_dir, model, tokenizer, encoders, model_name, args)
    print(f"\nFinal model ready for deployment at: {Path(save_dir).resolve()}")
    print("Point nlp_integration.py's MODEL_DIR at this folder.")


def build_arg_parser():
    p = argparse.ArgumentParser(description="Cross-validated benchmark and final trainer for the isolation-source classifier")
    p.add_argument("--data", required=True, help="Path to the labelled TSV")
    p.add_argument("--output", default="benchmark_results", help="Directory for comparison tables / per-class reports")
    p.add_argument("--models", nargs="+", default=list(models_to_benchmark.keys()),
                   choices=list(models_to_benchmark.keys()), help="Which encoders to benchmark")
    # Cross-validation and early stopping
    p.add_argument("--k_folds", type=int, default=5)
    p.add_argument("--max_epochs", type=int, default=15, help="Upper bound; early stopping usually stops sooner")
    p.add_argument("--patience", type=int, default=3, help="Stop after this many epochs without improvement")
    p.add_argument("--min_delta", type=float, default=0.001, help="Smallest score gain that counts as improvement")
    p.add_argument("--test_size", type=float, default=0.15, help="Fraction held out as the untouched test set")
    # Optimization
    p.add_argument("--batch_size", type=int, default=32)
    p.add_argument("--lr", type=float, default=2e-5)
    p.add_argument("--max_length", type=int, default=64)
    p.add_argument("--num_workers", type=int, default=2, help="DataLoader workers; raise on many-core CPUs")
    # Loss (long-tail handling)
    p.add_argument("--loss", choices=["ce", "focal"], default="ce",
                   help="ce = cross-entropy (default); focal focuses gradients on hard/rare classes")
    p.add_argument("--focal_gamma", type=float, default=2.0, help="Focal focusing strength; only used with --loss focal")
    p.add_argument("--no_class_weights", action="store_true",
                   help="Disable sqrt-inverse-frequency class weights (default: on)")
    # Final model
    p.add_argument("--train_final", default=None, choices=list(models_to_benchmark.keys()),
                   help="Skip benchmarking: train THIS model on ALL data and save")
    p.add_argument("--final_epochs", type=int, default=None,
                   help="Fixed epochs for --train_final (use the number the benchmark recommends)")
    p.add_argument("--save_dir", default="./models/isolation_source_v1",
                   help="Where --train_final saves the deployable model")
    # MLflow
    p.add_argument("--experiment_name", default="metaminer_isolation_source")
    p.add_argument("--mlflow_uri", default=None, help="Tracking URI; default is a local ./mlruns folder")
    p.add_argument("--no_mlflow", action="store_true", help="Disable MLflow logging entirely")
    return p


def main():
    args = build_arg_parser().parse_args()

    global mlflow_enabled
    mlflow_enabled = mlflow_available and not args.no_mlflow
    if mlflow_enabled:
        if args.mlflow_uri:
            mlflow.set_tracking_uri(args.mlflow_uri)
        mlflow.set_experiment(args.experiment_name)
        print(f"MLflow ON  (experiment='{args.experiment_name}', uri={mlflow.get_tracking_uri()})")
    elif not mlflow_available and not args.no_mlflow:
        print("MLflow not installed -- continuing without tracking (pip install mlflow to enable).")

    set_seed()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    if device.type == "cpu":
        print("NOTE: CPU training is slow; a GPU is strongly recommended.")
    print(f"Loss: {args.loss}"
          + (f" (gamma={args.focal_gamma})" if args.loss == "focal" else "")
          + f"  |  class weights: {'off' if args.no_class_weights else 'on'}")
    logging.info(f"Starting on {device} with loss={args.loss}, class_weights={not args.no_class_weights}.")

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = load_and_prepare(args.data)
    encoders = fit_label_encoders(df)
    encoder_lookup.update(encoders)

    # Final-training path: train the chosen model on everything and save.
    if args.train_final:
        with mlflow_run(run_name=f"final_{args.train_final}_{args.loss}"):
            train_final_model(args.train_final, df, encoders, args, device, args.save_dir)
        return

    # Benchmark path: hold out a test set, cross-validate the rest.
    dev_df, test_df = stratified_split(df, "Manual_identified_host", test_size=args.test_size)
    dev_df = dev_df.reset_index(drop=True)
    test_df = test_df.reset_index(drop=True)
    print(f"\nHeld-out split: dev={len(dev_df)}  test={len(test_df)} (test is touched only once, at the end)")

    all_aggs = {}
    for key in args.models:
        all_aggs[key] = run_cross_validation(key, models_to_benchmark[key], dev_df, encoders, args, device)

    if not all_aggs:
        print("No models completed.")
        return

    build_comparison(all_aggs, out_dir)
    winner = max(all_aggs, key=lambda k: all_aggs[k]["selection_score_mean"])
    rec_epochs = max(1, round(all_aggs[winner]["best_epoch_mean"]))
    print(f"\nRecommended model: {winner}  (CV selection score "
          f"{all_aggs[winner]['selection_score_mean']:.4f} +/- {all_aggs[winner]['selection_score_std']:.4f})")
    print(f"Recommended epochs for final training: {rec_epochs} (mean best epoch across folds)")
    logging.info(f"Recommended model: {winner} at ~{rec_epochs} epochs.")

    with mlflow_run(run_name=f"final_test_{winner}_{args.loss}"):
        mlflow_log_params({"winner": winner, "retrain_epochs": rec_epochs,
                           "n_dev": len(dev_df), "n_test": len(test_df)})
        retrain_and_test(winner, dev_df, test_df, encoders, args, device, rec_epochs, out_dir)

    print("\n" + "=" * 78)
    print("  NEXT STEP -- train the deployable model on ALL data:")
    print(f"    python {Path(__file__).name} --data {args.data} "
          f"--train_final {winner} --final_epochs {rec_epochs} --save_dir ./models/isolation_source_v1")
    print("=" * 78)
    print(f"\nAll tables/reports in {out_dir.resolve()}")
    if mlflow_enabled:
        print("Inspect curves and fold comparisons with:  mlflow ui")


if __name__ == "__main__":
    main()
