"""
Train a motif classifier using frozen ESM2-8M embeddings + XGBoost.
Usage: python training/train.py --motif lrr
Expects: training/data/lrr.fasta, training/data/lrr.labels
Outputs: models/lrr.ubj, models/lrr_meta.json
"""
import argparse
import json
import os
import numpy as np
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import classification_report, roc_auc_score, precision_recall_curve
from xgboost import XGBClassifier

from data import parse_fasta, parse_labels
from features import embed_sequences, make_windows

WINDOW_SIZE = 11
DATA_DIR    = "training/data"
MODELS_DIR  = "models"
ESM_MODEL   = "Synthyra/ESM2-8M"

MOTIF_CONFIG = {
    "LxxLxL"  : ("lrr.fasta", "lrr.labels",        6  + 10),
    "P-loop"  : ("nbs.fasta", "nbs_P-loop.labels",  9  + 10),
    "GLPL"    : ("nbs.fasta", "nbs_GLPL.labels",    5  + 10),
    "MHD"     : ("nbs.fasta", "nbs_MHD.labels",     3  + 10),
    "Walker-B": ("nbs.fasta", "nbs_Walker-B.labels", 8  + 10),
    "RNSB-A"  : ("nbs.fasta", "nbs_RNSB-A.labels",  10 + 10),
    "RNSB-B"  : ("nbs.fasta", "nbs_RNSB-B.labels",  7  + 10),
    "RNSB-C"  : ("nbs.fasta", "nbs_RNSB-C.labels",  10 + 10),
    "RNSB-D"  : ("nbs.fasta", "nbs_RNSB-D.labels",  9  + 10),
    "VG"      : ("nbs.fasta", "nbs_VG.labels",       5  + 10),
    "extEDVID": ("cc.fasta",  "cc_extEDVID.labels",  12 + 10),
    "aA"      : ("tir.fasta", "tir_aA.labels",       7  + 10),
    "aC"      : ("tir.fasta", "tir_aC.labels",       6  + 10),
    "aD3"     : ("tir.fasta", "tir_aD3.labels",      13 + 10),
    "bA"      : ("tir.fasta", "tir_bA.labels",       10 + 10),
    "bC"      : ("tir.fasta", "tir_bC.labels",       8  + 10),
    "bDaD1"   : ("tir.fasta", "tir_bDaD1.labels",    16 + 10),
}

def best_f1_threshold(y_true, y_proba):
    precisions, recalls, thresholds = precision_recall_curve(y_true, y_proba)
    f1_scores = 2 * (precisions * recalls) / (precisions + recalls + 1e-8)
    idx = np.argmax(f1_scores)
    return float(thresholds[idx]), float(precisions[idx]), float(recalls[idx]), float(f1_scores[idx])


def train(motif, fasta_path, labels_path, models_dir, window_size):
    # ── Load data ────────────────────────────────────────────────────────────
    print(f"Loading data for motif: {motif}")
    sequences_dict = parse_fasta(fasta_path)
    label_data     = parse_labels(labels_path)

    seq_ids   = [s for s in label_data if s in sequences_dict]
    sequences = [sequences_dict[s] for s in seq_ids]
    print(f"  {len(seq_ids)} sequences loaded")

    # ── Embed ────────────────────────────────────────────────────────────────
    print(f"Embedding with {ESM_MODEL}...")
    embeddings = embed_sequences(sequences, ESM_MODEL, device="cpu")

    # ── Build feature matrix ─────────────────────────────────────────────────
    print("Building feature matrix...")
    X_all, y_all, groups = [], [], []

    for seq_id, seq, emb in zip(seq_ids, sequences, embeddings):
        labels = label_data[seq_id]
        if emb.shape[0] != len(labels):
            print(f"  Warning: length mismatch for {seq_id} "
                  f"(emb={emb.shape[0]}, labels={len(labels)}), skipping")
            continue
        windows = make_windows(emb, window_size)
        y = np.array([l for _, _, l in labels], dtype=int)
        X_all.append(windows)
        y_all.append(y)
        groups.extend([seq_id] * len(y))

    X      = np.concatenate(X_all)
    y      = np.concatenate(y_all)
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
        n_estimators          = 500,
        max_depth             = 6,
        learning_rate         = 0.05,
        subsample             = 0.8,
        colsample_bytree      = 0.8,
        scale_pos_weight      = scale_pos_weight,
        eval_metric           = "logloss",
        early_stopping_rounds = 20,
        device                = "cpu",
        random_state          = 42,
    )
    model.fit(X_train, y_train, eval_set=[(X_val, y_val)], verbose=50)

    # ── Find best threshold ──────────────────────────────────────────────────
    y_proba = model.predict_proba(X_val)[:, 1]
    threshold, precision, recall, f1 = best_f1_threshold(y_val, y_proba)
    roc_auc = roc_auc_score(y_val, y_proba)

    print(f"\n── Validation results ──────────────────────────────────")
    print(f"  Best threshold : {threshold:.4f}")
    print(f"  Precision      : {precision:.4f}")
    print(f"  Recall         : {recall:.4f}")
    print(f"  F1             : {f1:.4f}")
    print(f"  ROC-AUC        : {roc_auc:.4f}")
    print()
    y_pred = (y_proba >= threshold).astype(int)
    print(classification_report(y_val, y_pred, target_names=["background", motif.upper()]))

    # ── Save model and metadata ──────────────────────────────────────────────
    os.makedirs(models_dir, exist_ok=True)
    model_path = os.path.join(models_dir, f"{motif}.ubj")
    meta_path  = os.path.join(models_dir, f"{motif}_meta.json")

    model.save_model(model_path)

    meta = {
        "motif"      : motif,
        "esm_model"  : ESM_MODEL,
        "window_size": window_size,
        "threshold"  : threshold,
        "val_metrics": {
            "precision": precision,
            "recall"   : recall,
            "f1"       : f1,
            "roc_auc"  : roc_auc,
        }
    }
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Model saved to {model_path}")
    print(f"Metadata saved to {meta_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--motif",  default=None, help="Single motif to train")
    parser.add_argument("--all",    action="store_true", help="Train all motifs")
    parser.add_argument("--fasta",  default=None, help="Override FASTA path")
    parser.add_argument("--labels", default=None, help="Override labels path")
    parser.add_argument("--data_dir",   default=DATA_DIR, help="Directory containing FASTA and label files")
    parser.add_argument("--models_dir", default=MODELS_DIR, help="Directory to save trained models")

    args = parser.parse_args()

    if args.all:
        for motif, (fasta, labels, window_size) in MOTIF_CONFIG.items():
            train(
                motif       = motif,
                fasta_path  = os.path.join(args.data_dir, fasta),
                labels_path = os.path.join(args.data_dir, labels),
                models_dir  = args.models_dir,
                window_size = window_size,
            )
    elif args.motif:
        fasta_path  = os.path.join(args.data_dir, MOTIF_CONFIG[args.motif][0])
        labels_path = os.path.join(args.data_dir, MOTIF_CONFIG[args.motif][1])
        window_size = MOTIF_CONFIG[args.motif][2]
        train(
            motif       = args.motif,
            fasta_path  = fasta_path,
            labels_path = labels_path,
            models_dir  = args.models_dir,
            window_size = window_size,
        )
    else:
        parser.error("Provide --motif or --all")