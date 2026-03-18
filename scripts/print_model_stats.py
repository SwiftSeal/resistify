#!/usr/bin/env python3
"""Print a markdown table of motif prediction model stats from model JSON files."""

import json
from pathlib import Path

models_dir = Path(__file__).parent.parent / "src/resistify/data/models"

rows = []
for meta_file in sorted(models_dir.glob("*_meta.json")):
    with open(meta_file) as f:
        data = json.load(f)
    motif = data["motif"]
    v = data["val_metrics"]
    rows.append((motif, v["precision"], v["recall"], v["f1"]))

header = "| Motif | Precision | Recall | F1 |"
sep = "| ----- | --------- | ------ | -- |"
print(header)
print(sep)
for motif, precision, recall, f1 in rows:
    print(f"| {motif} | {precision:.2f} | {recall:.2f} | {f1:.2f} |")
