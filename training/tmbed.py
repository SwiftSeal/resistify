"""Train a single TMBed CNN model using ESM2-8M embeddings.

Trains one Predictor model with an 80/20 train/val split (grouped by sequence),
using alpha-helix and signal-peptide data. Beta-barrel proteins are excluded.

Usage:
    uv run python training/tmbed.py [--device cpu|cuda] [--epochs 50]
                                    [--lr 1e-3] [--channels 64]
                                    [--output-dir src/resistify/data/tmbed_models]

Output: <output-dir>/tmbed_model.pt
"""

import argparse
import logging
import random
import sys
from pathlib import Path

import torch
import torch.nn.functional as F
from sklearn.model_selection import GroupShuffleSplit
from torch.optim import Adam
from transformers import AutoModel
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))
from resistify.tmbed import Predictor

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "TMbed" / "data" / "datasets"
OUTPUT_DIR = Path(__file__).parent.parent / "src" / "resistify" / "data" / "tmbed_models"

ESM_MODEL = "Synthyra/ESM2-8M"
ESM_REVISION = "f3c6441"

LABEL_MAP = {
    "H": 0,
    "h": 0,
    "S": 1,
    "1": 2,
    "2": 3,
}
CLASS_NAMES = ["H/h (TM helix)", "S (signal peptide)", "1 (inside)", "2 (outside)"]
IGNORE_INDEX = -100
BETA_CHARS = set("Bb")


def parse_fasta(path: Path):
    """Parse 3-line FASTA (header / sequence / annotation). Returns list of (header, seq, ann)."""
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


def filter_beta(records):
    kept = [(h, s, a) for h, s, a in records if not BETA_CHARS.intersection(a)]
    skipped = len(records) - len(kept)
    if skipped:
        logger.info(f"  Skipped {skipped} beta-barrel sequences")
    return kept


def ann_to_labels(ann: str) -> torch.Tensor:
    labels = torch.full((len(ann),), IGNORE_INDEX, dtype=torch.long)
    for i, c in enumerate(ann):
        if c in LABEL_MAP:
            labels[i] = LABEL_MAP[c]
    return labels


def load_esm(device: str):
    logger.info(f"Loading {ESM_MODEL}...")
    model = (
        AutoModel.from_pretrained(
            ESM_MODEL,
            trust_remote_code=True,
            revision=ESM_REVISION,
        )
        .eval()
        .to(device)
    )
    logger.info("ESM2-8M loaded.")
    return model, model.tokenizer


@torch.inference_mode()
def embed_sequence(esm_model, tokenizer, seq: str, device: str) -> torch.Tensor:
    """Returns float32 tensor of shape (L, 320)."""
    encoded = tokenizer(seq, return_tensors="pt").to(device)
    out = esm_model(**encoded)
    return out.last_hidden_state[0, 1:-1, :].float().cpu()


def embed_all(records, esm_model, tokenizer, device: str):
    """Pre-compute embeddings. Returns list of (header, emb, labels)."""
    dataset = []
    for header, seq, ann in tqdm(records, desc="  Embedding", leave=False):
        emb = embed_sequence(esm_model, tokenizer, seq, device)
        labels = ann_to_labels(ann)
        assert emb.shape[0] == len(seq), f"Embedding length mismatch for {header}"
        dataset.append((header, emb, labels))
    return dataset


def compute_class_weights(dataset, num_classes: int = 4) -> torch.Tensor:
    """Inverse-frequency class weights computed from training labels."""
    counts = torch.zeros(num_classes, dtype=torch.long)
    for _, _, labels in dataset:
        valid = labels[labels != IGNORE_INDEX]
        for cls in range(num_classes):
            counts[cls] += (valid == cls).sum()
    total = counts.sum().float()
    weights = total / (num_classes * counts.float())
    return weights


def train_epoch(model, dataset, optimizer, device: str, class_weights: torch.Tensor):
    model.train()
    total_loss = 0.0
    n = 0
    w = class_weights.to(device)
    order = list(range(len(dataset)))
    random.shuffle(order)
    for i in order:
        _, emb, labels = dataset[i]
        L = emb.shape[0]
        x = emb.unsqueeze(0).to(device)
        mask = torch.ones(1, L, dtype=x.dtype, device=device)
        y = labels.to(device)

        logits = model(x, mask).squeeze(0).permute(1, 0)  # (L, 4)
        loss = F.cross_entropy(logits, y, weight=w, ignore_index=IGNORE_INDEX)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        total_loss += loss.item()
        n += 1
    return total_loss / n if n > 0 else 0.0


@torch.no_grad()
def eval_epoch(model, dataset, device: str, class_weights: torch.Tensor):
    model.eval()
    total_loss = 0.0
    n = 0
    w = class_weights.to(device)
    for _, emb, labels in dataset:
        L = emb.shape[0]
        x = emb.unsqueeze(0).to(device)
        mask = torch.ones(1, L, dtype=x.dtype, device=device)
        y = labels.to(device)

        logits = model(x, mask).squeeze(0).permute(1, 0)  # (L, 4)
        loss = F.cross_entropy(logits, y, weight=w, ignore_index=IGNORE_INDEX)

        total_loss += loss.item()
        n += 1
    return total_loss / n if n > 0 else 0.0


@torch.no_grad()
def compute_accuracy(model, dataset, device: str):
    """Compute overall and per-class per-residue accuracy on a dataset."""
    model.eval()

    num_classes = 4
    correct = torch.zeros(num_classes, dtype=torch.long)
    total = torch.zeros(num_classes, dtype=torch.long)
    overall_correct = 0
    overall_total = 0

    for _, emb, labels in dataset:
        L = emb.shape[0]
        x = emb.unsqueeze(0).to(device)
        mask = torch.ones(1, L, dtype=x.dtype, device=device)
        y = labels  # (L,) on CPU

        logits = model(x, mask).squeeze(0).permute(1, 0).cpu()  # (L, 4)
        preds = logits.argmax(dim=1)  # (L,)

        valid = y != IGNORE_INDEX
        y_valid = y[valid]
        preds_valid = preds[valid]

        overall_correct += (preds_valid == y_valid).sum().item()
        overall_total += valid.sum().item()

        for cls in range(num_classes):
            mask_cls = y_valid == cls
            if mask_cls.any():
                correct[cls] += (preds_valid[mask_cls] == cls).sum()
                total[cls] += mask_cls.sum()

    return overall_correct, overall_total, correct, total


def main():
    parser = argparse.ArgumentParser(description="Train single TMBed CNN with ESM2-8M")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--channels", type=int, default=64)
    parser.add_argument("--patience", type=int, default=10)
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))
    args = parser.parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Device: {args.device}, Epochs: {args.epochs}, LR: {args.lr}, Channels: {args.channels}")

    # ── Load and merge datasets ───────────────────────────────────────────────
    logger.info("Loading datasets...")
    alpha_records = filter_beta(parse_fasta(DATA_DIR / "alpha.fasta"))
    signalp_records = filter_beta(parse_fasta(DATA_DIR / "signalp.fasta"))
    logger.info(f"  alpha.fasta:   {len(alpha_records)} sequences")
    logger.info(f"  signalp.fasta: {len(signalp_records)} sequences")

    records = alpha_records + signalp_records
    # Tag each record with its source so we can oversample after embedding
    sources = ["alpha"] * len(alpha_records) + ["signalp"] * len(signalp_records)

    # ── Train / val split grouped by sequence ─────────────────────────────────
    headers = [h for h, _, _ in records]
    splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    train_idx, val_idx = next(splitter.split(records, groups=headers))

    train_records = [records[i] for i in train_idx]
    train_sources = [sources[i] for i in train_idx]
    val_records = [records[i] for i in val_idx]
    logger.info(f"  Train: {len(train_records)} sequences, Val: {len(val_records)} sequences")

    # ── Embed ─────────────────────────────────────────────────────────────────
    esm_model, tokenizer = load_esm(args.device)

    logger.info("Embedding training sequences...")
    train_data = embed_all(train_records, esm_model, tokenizer, args.device)

    logger.info("Embedding validation sequences...")
    val_data = embed_all(val_records, esm_model, tokenizer, args.device)

    # ── Oversample alpha sequences to match signalp count ─────────────────────
    alpha_data = [s for s, src in zip(train_data, train_sources) if src == "alpha"]
    signalp_data = [s for s, src in zip(train_data, train_sources) if src == "signalp"]
    if alpha_data and signalp_data:
        repeat = max(1, round(len(signalp_data) / len(alpha_data)))
        train_data = signalp_data + alpha_data * repeat
        logger.info(f"  Oversampled alpha {repeat}x: {len(alpha_data)} → {len(alpha_data) * repeat} sequences")
        logger.info(f"  Effective train size: {len(train_data)} sequences")

    # ── Class weights (inverse frequency from training data) ──────────────────
    class_weights = compute_class_weights(train_data)
    logger.info("  Class weights: " + ", ".join(
        f"{CLASS_NAMES[i].split()[0]}={class_weights[i]:.2f}" for i in range(4)
    ))

    # ── Train ─────────────────────────────────────────────────────────────────
    model = Predictor(channels=args.channels).to(args.device)
    optimizer = Adam(model.parameters(), lr=args.lr)

    best_val_loss = float("inf")
    best_state = None
    patience_count = 0

    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(model, train_data, optimizer, args.device, class_weights)
        val_loss = eval_epoch(model, val_data, args.device, class_weights)

        if epoch % 5 == 0 or epoch == 1:
            logger.info(f"  Epoch {epoch:3d}/{args.epochs}  train={train_loss:.4f}  val={val_loss:.4f}")

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_count = 0
        else:
            patience_count += 1
            if patience_count >= args.patience:
                logger.info(f"  Early stopping at epoch {epoch} (patience={args.patience})")
                break

    logger.info(f"  Best val loss: {best_val_loss:.4f}")

    # ── Validation accuracy ───────────────────────────────────────────────────
    model.load_state_dict(best_state)
    overall_correct, overall_total, correct, total = compute_accuracy(model, val_data, args.device)

    logger.info("\n── Validation accuracy ──────────────────────────────────")
    logger.info(f"  Overall: {overall_correct}/{overall_total} = {overall_correct/overall_total:.4f}")
    for cls in range(4):
        if total[cls] > 0:
            acc = correct[cls].item() / total[cls].item()
            logger.info(f"  {CLASS_NAMES[cls]:<22}: {correct[cls]:6d}/{total[cls]:6d} = {acc:.4f}")
        else:
            logger.info(f"  {CLASS_NAMES[cls]:<22}: no examples")

    # ── Save ──────────────────────────────────────────────────────────────────
    out_path = Path(args.output_dir) / "tmbed_model.pt"
    torch.save({"model": best_state}, out_path)
    logger.info(f"\nModel saved to {out_path}")


if __name__ == "__main__":
    main()
