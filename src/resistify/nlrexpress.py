import json
import logging
import numpy as np
import torch
from pathlib import Path
from transformers import AutoModel, AutoTokenizer
from xgboost import XGBClassifier
from resistify.annotation import Protein, Annotation

logger = logging.getLogger(__name__)

ESM_MODEL = "Synthyra/ESM2-8M"
MODELS_DIR = Path(__file__).parent / "data" / "models"

WINDOW_SIZE = 30

MOTIF_SPAN_LENGTHS = {
    "VG": 5,
    "P-loop": 9,
    "RNBS-A": 10,
    "RNBS-B": 7,
    "RNBS-C": 10,
    "RNBS-D": 9,
    "Walker-B": 8,
    "GLPL": 5,
    "MHD": 3,
    "extEDVID": 12,
    "aA": 7,
    "aC": 6,
    "aD3": 13,
    "bA": 10,
    "bC": 8,
    "bDaD1": 16,
    "LxxLxL": 6,
}


# ── Model loading ─────────────────────────────────────────────────────────────


def _load_models(models_dir: Path, search_type: str, threads: int):
    models = {}
    motifs = list(MOTIF_SPAN_LENGTHS.keys()) if search_type == "all" else [search_type]

    for motif in motifs:
        model_path = models_dir / f"{motif}.ubj"
        meta_path = models_dir / f"{motif}_meta.json"

        if not model_path.exists():
            logger.warning(f"No model found for motif {motif}, skipping")
            continue

        clf = XGBClassifier()
        clf.load_model(model_path)
        clf.set_params(nthread=threads)

        with open(meta_path) as f:
            meta = json.load(f)

        models[motif] = {"model": clf, "meta": meta}
        logger.debug(f"Loaded {motif} (threshold={meta['threshold']:.4f})")

    return models


# ── ESM embedding ─────────────────────────────────────────────────────────────


def _load_esm(device: str) -> tuple[AutoModel, AutoTokenizer]:
    logger.info(f"Loading {ESM_MODEL} on {device}...")
    tokenizer = AutoTokenizer.from_pretrained(ESM_MODEL, trust_remote_code=True)
    model = (
        AutoModel.from_pretrained(
            ESM_MODEL,
            trust_remote_code=True,
        )
        .eval()
        .to(device)
    )
    return model, tokenizer


def _embed(model, tokenizer, sequence: str, device: str) -> np.ndarray:
    inputs = tokenizer(sequence, return_tensors="pt").to(device)
    with torch.no_grad():
        outputs = model(**inputs)
    return outputs.last_hidden_state[0, 1:-1].cpu().float().numpy()


# ── Feature extraction ────────────────────────────────────────────────────────


def _make_windows(matrix: np.ndarray, window_size: int) -> np.ndarray:
    pad = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    return np.stack([padded[i : i + window_size].flatten() for i in range(len(matrix))])


# ── Main entry point ──────────────────────────────────────────────────────────


def nlrexpress(
    proteins: dict[str, Protein],
    search_type: str = "all",
    device: str = "cpu",
    threads: int = 1,
):
    logger.info(f"Running motif classifier for '{search_type}' motifs")

    models = _load_models(MODELS_DIR, search_type, threads)
    if not models:
        logger.warning("No models loaded, skipping")
        return proteins

    torch.set_num_threads(threads)
    esm, tokenizer = _load_esm(device)

    for seq_id, protein in proteins.items():
        emb = _embed(esm, tokenizer, protein.sequence, device)

        if len(emb) < WINDOW_SIZE:
            continue

        windows = _make_windows(emb, WINDOW_SIZE)

        for motif, clf in models.items():
            threshold = clf["meta"]["threshold"]
            span = MOTIF_SPAN_LENGTHS[motif]

            proba = clf["model"].predict_proba(windows)[:, 1]

            for idx in np.where(proba >= threshold)[0]:
                end = int(idx + span)
                if end > protein.length:
                    continue
                protein.add_annotation(
                    Annotation(
                        name=motif,
                        type="motif",
                        start=int(idx + 1),
                        end=end,
                        source="motif_classifier",
                        score=float(proba[idx]),
                    )
                )

    logger.info("Motif classification completed")
    return proteins
