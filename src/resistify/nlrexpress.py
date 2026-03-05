import json
import logging
import sqlite3
import tempfile
import numpy as np
import torch
from pathlib import Path
from tqdm.auto import tqdm
from transformers import AutoModel
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


def _load_esm(device: str) -> AutoModel:
    logger.info(f"Loading {ESM_MODEL} on {device}...")
    model = (
        AutoModel.from_pretrained(
            ESM_MODEL,
            trust_remote_code=True,
        )
        .eval()
        .to(device)
    )
    return model


def _embed_all(model, sequences: list[str], db_path: str):
    model.embed_dataset(
        sequences=sequences,
        tokenizer=model.tokenizer,
        batch_size=1,  # Batching is broken but 1 is safe - revisit
        max_len=None,
        full_embeddings=True,
        embed_dtype=torch.float32,
        num_workers=0,
        sql=True,
        sql_db_path=db_path,
        save=False,
    )


def _fetch_embedding(conn: sqlite3.Connection, sequence: str) -> np.ndarray | None:
    row = conn.execute(
        "SELECT embedding, shape, dtype FROM embeddings WHERE sequence = ?", (sequence,)
    ).fetchone()
    if row is None:
        return None
    blob, shape_str, dtype_str = row
    shape = tuple(int(x) for x in shape_str.strip("()").split(",") if x.strip())
    emb = np.frombuffer(blob, dtype=np.dtype(dtype_str)).reshape(shape).copy()
    if emb.shape[0] == len(sequence) + 2:
        emb = emb[1:-1]
    if np.isnan(emb).any():
        raise RuntimeError(f"Empty embeddings found for sequence {sequence}")
    return emb


def _make_windows(matrix: np.ndarray, window_size: int) -> np.ndarray:
    pad = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    windowed = np.lib.stride_tricks.sliding_window_view(padded, window_size, axis=0)[
        : len(matrix)
    ]
    return windowed.transpose(0, 2, 1).reshape(len(matrix), -1)


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
    esm = _load_esm(device)

    sequences = [p.sequence for p in proteins.values() if p.length >= WINDOW_SIZE]
    logger.info(f"Embedding {len(sequences)} sequences...")

    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = f.name
    try:
        _embed_all(esm, sequences, db_path)
        conn = sqlite3.connect(db_path)

        for seq_id, protein in tqdm(proteins.items(), desc="Predicting motifs"):
            if protein.length < WINDOW_SIZE:
                continue

            emb = _fetch_embedding(conn, protein.sequence)

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

        conn.close()
    finally:
        Path(db_path).unlink(missing_ok=True)

    logger.info("Motif classification completed")
    return proteins
