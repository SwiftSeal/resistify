"""
Train a motif classifier using frozen ESM2-8M embeddings + XGBoost.
Usage: python training/train.py --motif extEDVID
Expects: training/data/<motif>.labels
Outputs: src/resistify/data/models/<motif>.ubj, <motif>_meta.json
"""

import gc
import os
import argparse
import json
import numpy as np
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import classification_report, roc_auc_score, precision_recall_curve
from xgboost import XGBClassifier

from data import parse_labels
from features import embed_sequences, load_esm_model, make_windows

WINDOW_SIZE = 11
DATA_DIR = "training/data"
MODELS_DIR = "src/resistify/data/models"
ESM_MODEL = "Synthyra/ESM2-8M"

MOTIF_CONFIG = {
    "extEDVID": ("extEDVID.labels", 30),
    "P-loop":   ("P-loop.labels",   30),
    "GLPL":     ("GLPL.labels",     30),
    "MHD":      ("MHD.labels",      30),
    "Walker-B": ("Walker-B.labels", 30),
    "RNBS-A":   ("RNBS-A.labels",  30),
    "RNBS-B":   ("RNBS-B.labels",   30),
    "RNBS-C":   ("RNBS-C.labels",  30),
    "RNBS-D":   ("RNBS-D.labels",   30),
    "VG":       ("VG.labels",       30),
    # "aA":    ("aA.labels",    7 + 10),
    # "aC":    ("aC.labels",    6 + 10),
    # "aD3":   ("aD3.labels",  13 + 10),
    # "bA":    ("bA.labels",   10 + 10),
    # "bC":    ("bC.labels",    8 + 10),
    # "bDaD1": ("bDaD1.labels",16 + 10),
    "LxxLxL":   ("LxxLxL.labels",   30),
}


def best_f1_threshold(y_true, y_proba):
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_proba)
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)
    idx = np.argmax(f1_scores)
    return (
        min(float(thresholds[idx]), 0.99),
        float(precisions[idx]),
        float(recalls[idx]),
        float(f1_scores[idx]),
    )


def train(motif, labels_path, models_dir, window_size, esm_model=None):
    # ── Load data ────────────────────────────────────────────────────────────
    print(f"Loading data for motif: {motif}")
    label_data = parse_labels(labels_path)

    seq_ids = list(label_data.keys())
    sequences = ["".join(r for _, r, _ in label_data[s]) for s in seq_ids]
    print(f"  {len(seq_ids)} sequences loaded")

    # ── Embed ────────────────────────────────────────────────────────────────
    print(f"Embedding with {ESM_MODEL}...")
    embeddings = embed_sequences(sequences, ESM_MODEL, device="cpu", model=esm_model)

    # ── Build feature matrix ─────────────────────────────────────────────────
    print("Building feature matrix...")
    X_all, y_all, groups = [], [], []

    for seq_id, seq, emb in zip(seq_ids, sequences, embeddings):
        labels = label_data[seq_id]
        if emb.shape[0] != len(labels):
            print(
                f"  Warning: length mismatch for {seq_id} "
                f"(emb={emb.shape[0]}, labels={len(labels)}), skipping"
            )
            continue
        windows = make_windows(emb, window_size)
        y = np.array([label for _, _, label in labels], dtype=int)
        X_all.append(windows)
        y_all.append(y)
        groups.extend([seq_id] * len(y))

    X = np.concatenate(X_all)
    y = np.concatenate(y_all)
    groups = np.array(groups)

    n_pos = y.sum()
    n_neg = (y == 0).sum()
    print(f"  Total positions: {len(y)}  |  Positive: {n_pos}  |  Negative: {n_neg}")

    # ── Train / val split ────────────────────────────────────────────────────
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    train_idx, val_idx = next(splitter.split(X, y, groups=groups))
    X_train, X_val = X[train_idx], X[val_idx]
    y_train, y_val = y[train_idx], y[val_idx]

    scale_pos_weight = (y_train == 0).sum() / (y_train == 1).sum()
    print(f"  scale_pos_weight: {scale_pos_weight:.1f}")

    # ── Fit XGBoost ──────────────────────────────────────────────────────────
    print("Training XGBoost...")
    model = XGBClassifier(
        n_estimators=500,
        max_depth=6,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=scale_pos_weight,
        eval_metric="logloss",
        early_stopping_rounds=20,
        device="cpu",
        random_state=42,
    )
    model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=50)

    # ── Find best threshold ──────────────────────────────────────────────────
    y_proba = model.predict_proba(X_val)[:, 1]
    threshold, precision, recall, f1 = best_f1_threshold(y_val, y_proba)
    roc_auc = roc_auc_score(y_val, y_proba)

    print("\n── Validation results ──────────────────────────────────")
    print(f"  Best threshold : {threshold:.4f}")
    print(f"  Precision      : {precision:.4f}")
    print(f"  Recall         : {recall:.4f}")
    print(f"  F1             : {f1:.4f}")
    print(f"  ROC-AUC        : {roc_auc:.4f}")
    print()
    y_pred = (y_proba >= threshold).astype(int)
    print(
        classification_report(y_val, y_pred, target_names=["background", motif.upper()])
    )

    # ── Save model and metadata ──────────────────────────────────────────────
    os.makedirs(models_dir, exist_ok=True)
    model_path = os.path.join(models_dir, f"{motif}.ubj")
    meta_path = os.path.join(models_dir, f"{motif}_meta.json")

    model.save_model(model_path)

    meta = {
        "motif": motif,
        "esm_model": ESM_MODEL,
        "window_size": window_size,
        "threshold": threshold,
        "val_metrics": {
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "roc_auc": roc_auc,
        },
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Model saved to {model_path}")
    print(f"Metadata saved to {meta_path}")

    return {
        "motif": motif,
        "threshold": threshold,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "roc_auc": roc_auc,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--motif", default=None, help="Single motif to train")
    parser.add_argument("--all", action="store_true", help="Train all motifs")
    parser.add_argument("--labels", default=None, help="Override labels path")
    parser.add_argument(
        "--data_dir",
        default=DATA_DIR,
        help="Directory containing label files",
    )
    parser.add_argument(
        "--models_dir", default=MODELS_DIR, help="Directory to save trained models"
    )

    args = parser.parse_args()

    if args.all:
        esm_model = load_esm_model(ESM_MODEL, device="cpu")
        results = []
        for motif, (labels, window_size) in MOTIF_CONFIG.items():
            result = train(
                motif=motif,
                labels_path=os.path.join(args.data_dir, labels),
                models_dir=args.models_dir,
                window_size=window_size,
                esm_model=esm_model,
            )
            results.append(result)
            gc.collect()

        print("\n── All-motif summary ───────────────────────────────────")
        print(f"{'Motif':<12}  {'Threshold':>9}  {'Precision':>9}  {'Recall':>9}  {'F1':>9}  {'ROC-AUC':>9}")
        print("-" * 67)
        for r in results:
            print(
                f"{r['motif']:<12}  {r['threshold']:>9.4f}  {r['precision']:>9.4f}"
                f"  {r['recall']:>9.4f}  {r['f1']:>9.4f}  {r['roc_auc']:>9.4f}"
            )
    elif args.motif:
        labels_path = args.labels or os.path.join(args.data_dir, MOTIF_CONFIG[args.motif][0])
        window_size = MOTIF_CONFIG[args.motif][1]
        train(
            motif=args.motif,
            labels_path=labels_path,
            models_dir=args.models_dir,
            window_size=window_size,
        )
    else:
        parser.error("Provide --motif or --all")
