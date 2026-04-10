"""Evaluate the trained TMBed model against the validation split.

Recreates the same 80/20 train/val split used during training, runs predictions
through the full inference pipeline (ESM2 → CNN → CRF decoder), and writes a
comparison file showing original vs predicted annotations per sequence.

Usage:
    uv run python training/evaluate_tmbed.py [--device cpu|cuda] [--output predictions.txt]

Output format (one block per sequence):
    >header
    Seq:  MYGKIIFVLL...
    Orig: SSSSSSHHHH...   (original annotation from dataset)
    Pred: SSSSSShhhh...   (model prediction via CRF decoder)
"""

import argparse
import sys
from pathlib import Path

import torch
from sklearn.model_selection import GroupShuffleSplit
from transformers import AutoModel

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from resistify.tmbed import Decoder, Predictor, PRED_MAP

DATA_DIR = Path(__file__).parent.parent / "TMbed" / "data" / "datasets"
MODEL_PATH = Path(__file__).parent.parent / "src" / "resistify" / "data" / "tmbed_models" / "tmbed_model.pt"

ESM_MODEL = "Synthyra/ESM2-8M"
ESM_REVISION = "f3c6441"

BETA_CHARS = set("Bb")


def parse_fasta(path: Path):
    records = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if not lines[i].startswith(">"):
            i += 1
            continue
        if i + 2 >= len(lines):
            break
        header = lines[i][1:].strip()
        seq = lines[i + 1].strip()
        ann = lines[i + 2].strip()
        if ann.startswith(">") or len(ann) != len(seq):
            i += 2
            continue
        records.append((header, seq, ann))
        i += 3
    return records


def load_val_records():
    alpha = [r for r in parse_fasta(DATA_DIR / "alpha.fasta") if not BETA_CHARS.intersection(r[2])]
    signalp = [r for r in parse_fasta(DATA_DIR / "signalp.fasta") if not BETA_CHARS.intersection(r[2])]
    records = alpha + signalp

    headers = [h for h, _, _ in records]
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    _, val_idx = next(splitter.split(records, groups=headers))
    return [records[i] for i in val_idx]


@torch.inference_mode()
def embed_sequence(esm_model, tokenizer, seq: str, device: str) -> tuple[torch.Tensor, torch.Tensor]:
    """Returns (embedding, attention_mask) with BOS stripped, shape (1, L, C) and (1, L)."""
    aa_map = str.maketrans("BJOUZ", "XXXXX")
    seq_clean = seq.upper().translate(aa_map)
    encoded = tokenizer(seq_clean, return_tensors="pt").to(device)
    out = esm_model(**encoded)
    emb = out.last_hidden_state[:, 1:, :].float()
    mask = encoded["attention_mask"][:, 1:].float()
    return emb, mask


def predict(esm_model, tokenizer, predictor, decoder, seq: str, device: str) -> str:
    emb, mask = embed_sequence(esm_model, tokenizer, seq, device)
    with torch.no_grad():
        logits = predictor(emb, mask)
        probs = torch.softmax(logits, dim=1)
    decoded = decoder(probs.cpu(), mask.cpu()).byte()
    seq_len = len(seq)
    return "".join(PRED_MAP[int(x)] for x in decoded[0, :seq_len])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--output", default="predictions.txt")
    args = parser.parse_args()

    print("Loading validation split...")
    val_records = load_val_records()
    print(f"  {len(val_records)} sequences")

    print(f"Loading ESM2 ({ESM_MODEL})...")
    esm_model = (
        AutoModel.from_pretrained(ESM_MODEL, trust_remote_code=True, revision=ESM_REVISION)
        .eval()
        .to(args.device)
    )
    tokenizer = esm_model.tokenizer

    print(f"Loading predictor from {MODEL_PATH}...")
    predictor = Predictor()
    predictor.load_state_dict(torch.load(MODEL_PATH, weights_only=True)["model"])
    predictor = predictor.eval().to(args.device)

    decoder = Decoder()

    out_path = Path(args.output)
    print(f"Writing predictions to {out_path}...")

    with out_path.open("w") as f:
        for i, (header, seq, orig_ann) in enumerate(val_records, 1):
            pred_ann = predict(esm_model, tokenizer, predictor, decoder, seq, args.device)
            f.write(f">{header}\n")
            f.write(f"Seq:  {seq}\n")
            f.write(f"Orig: {orig_ann}\n")
            f.write(f"Pred: {pred_ann}\n")
            f.write("\n")

            if i % 100 == 0:
                print(f"  {i}/{len(val_records)}")

    print(f"Done. {len(val_records)} sequences written to {out_path}")


if __name__ == "__main__":
    main()
