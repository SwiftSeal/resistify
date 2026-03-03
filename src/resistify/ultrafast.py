import json
import logging
import sqlite3
import pickle
import numpy as np
import torch
from pathlib import Path
from transformers import AutoModel
from xgboost import XGBClassifier
from resistify.annotation import Protein, Annotation

logger = logging.getLogger(__name__)

ESM_MODEL = "Synthyra/ESM2-8M"

MOTIF_SPAN_LENGTHS = {
    "LxxLxL" : 6,
    "VG"     : 5,
    "P-loop" : 9,
    "RNSB-A" : 10,
    "RNSB-B" : 7,
    "RNSB-C" : 10,
    "RNSB-D" : 9,
    "Walker-B": 8,
    "GLPL"   : 5,
    "MHD"    : 3,
    "extEDVID": 12,
    "aA"     : 7,
    "aC"     : 6,
    "aD3"    : 13,
    "bA"     : 10,
    "bC"     : 8,
    "bDaD1"  : 16,
}


# ── Model loading ─────────────────────────────────────────────────────────────

def _load_models(models_dir: Path, search_type: str):
    models = {}
    motifs = list(MOTIF_SPAN_LENGTHS.keys()) if search_type == "all" else [search_type]

    for motif in motifs:
        model_path = models_dir / f"{motif}.ubj"
        meta_path  = models_dir / f"{motif}_meta.json"

        if not model_path.exists():
            logger.warning(f"No model found for motif {motif}, skipping")
            continue

        clf = XGBClassifier()
        clf.load_model(model_path)

        with open(meta_path) as f:
            meta = json.load(f)

        models[motif] = {"model": clf, "meta": meta}
        logger.debug(f"Loaded {motif} (threshold={meta['threshold']:.4f})")

    return models


# ── ESM embedding ─────────────────────────────────────────────────────────────

def _load_esm(device):
    logger.info(f"Loading {ESM_MODEL} on {device}...")
    model = AutoModel.from_pretrained(
        ESM_MODEL,
        trust_remote_code=True,
    ).half().eval().to(device)
    return model


def _embed_and_save(sequences: list[str], db_path: Path, device: str):
    """Embed all sequences and store in SQLite."""
    esm = _load_esm(device)
    logger.info(f"Embedding {len(sequences)} sequences -> {db_path}")
    esm.embed_dataset(
        sequences       = sequences,
        tokenizer       = esm.tokenizer,
        batch_size      = 16,
        max_len         = None,
        full_embeddings = True,
        embed_dtype     = torch.float32,
        num_workers     = 0,
        sql             = True,
        sql_db_path     = str(db_path),
        save            = False,
    )


def _load_embedding(conn: sqlite3.Connection, sequence: str) -> np.ndarray | None:
    """Retrieve a single embedding from SQLite by sequence string."""
    row = conn.execute(
        "SELECT embedding, shape, dtype FROM embeddings WHERE sequence = ?",
        (sequence,)
    ).fetchone()

    if row is None:
        return None

    blob, shape_str, dtype_str = row
    shape = tuple(int(x) for x in shape_str.strip("()").split(",") if x.strip())
    emb   = np.frombuffer(blob, dtype=np.dtype(dtype_str)).reshape(shape).copy()

    # Strip BOS/EOS tokens if present
    if emb.shape[0] == len(sequence) + 2:
        emb = emb[1:-1]

    return emb


# ── Feature extraction ────────────────────────────────────────────────────────

def _make_windows(matrix: np.ndarray, window_size: int) -> np.ndarray:
    pad    = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    return np.stack([padded[i:i + window_size].flatten() for i in range(len(matrix))])


# ── Main entry point ──────────────────────────────────────────────────────────

def motif_classifier(
    proteins: dict[str, Protein],
    models_dir: Path,
    search_type: str = "all",
    device: str = "cpu",
    db_path: Path | None = None,
):
    """
    Run motif classification on a dict of Protein objects.

    Args:
        proteins:    dict of {seq_id: Protein}
        models_dir:  directory containing .ubj and _meta.json model files
        search_type: "all" or a specific motif name
        device:      "cpu" or "cuda"
        db_path:     optional path to pre-computed SQLite embeddings.
                     If None, embeddings are computed on the fly and discarded.
    """
    logger.info(f"Running motif classifier for '{search_type}' motifs")

    models = _load_models(models_dir, search_type)
    if not models:
        logger.warning("No models loaded, skipping")
        return proteins

    sequences    = [p.sequence for p in proteins.values()]
    sequence_ids = list(proteins.keys())

    # Embed if no db provided, otherwise verify db exists
    if db_path is None:
        import tempfile, atexit, os
        tmp = tempfile.NamedTemporaryFile(suffix=".db", delete=False)
        tmp.close()
        db_path = Path(tmp.name)
        atexit.register(os.unlink, db_path)
        _embed_and_save(sequences, db_path, device)
    elif not db_path.exists():
        logger.info(f"No embedding DB found at {db_path}, computing...")
        _embed_and_save(sequences, db_path, device)
    else:
        logger.info(f"Using existing embeddings from {db_path}")

    # Classify per sequence
    conn = sqlite3.connect(db_path)

    for seq_id, seq in zip(sequence_ids, sequences):
        emb = _load_embedding(conn, seq)

        if emb is None:
            logger.warning(f"No embedding found for {seq_id}, skipping")
            continue

        for motif, clf in models.items():
            window_size = clf["meta"]["window_size"]
            threshold   = clf["meta"]["threshold"]
            span        = MOTIF_SPAN_LENGTHS[motif]

            if len(emb) < window_size:
                continue

            windows = _make_windows(emb, window_size)
            proba   = clf["model"].predict_proba(windows)[:, 1]

            for idx in np.where(proba >= threshold)[0]:
                proteins[seq_id].add_annotation(
                    Annotation(
                        name   = motif,
                        type   = "motif",
                        start  = int(idx + 1),
                        end    = int(idx + span),
                        source = "motif_classifier",
                        score  = float(proba[idx]),
                    )
                )

    conn.close()
    logger.info("Motif classification completed")
    return proteins