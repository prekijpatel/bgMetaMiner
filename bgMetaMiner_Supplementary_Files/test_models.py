from isolation_source_model import IsolationSourceModel

model = IsolationSourceModel.load("./models/isolation_source_v1")

examples = [
    ("Homo sapiens", "UTI", "urine sample"),
    ("Infant", "Severe sepsis", "blood"),
    (None, None, "ventilator tube tip"),
    ("Dog", None, "wound swab"),
    (None, None, "municipal wastewater sludge"),
    ("Synthetic variant of Salmonella Typhi", None, None),
]

# predict_batch takes the whole list and returns [(labels, confs), ...]
for (host, disease, source), (labels, confs) in zip(examples,
        model.predict_batch(examples, return_confidence=True)):
    print(f"\nInput : host={host!r}  disease={disease!r}  source={source!r}")
    print(f"  host     = {labels[0]:35s} (conf {confs[0]:.2f})")
    print(f"  category = {labels[1]:35s} (conf {confs[1]:.2f})")
    print(f"  source   = {labels[2]:35s} (conf {confs[2]:.2f})")
    print(f"  sample   = {labels[3]:35s} (conf {confs[3]:.2f})")