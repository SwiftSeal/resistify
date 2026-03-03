"""
Feature extraction utilities for motif classifier training.
"""

import sqlite3
import tempfile
import numpy as np
import torch
from pathlib import Path
from transformers import AutoModel


def load_esm_model(model_path, device):
    """Load a frozen ESM2 model onto the specified device."""
    print(f"Loading {model_path} on {device}...")
    model = (
        AutoModel.from_pretrained(
            model_path,
            trust_remote_code=True,
        )
        .eval()
        .to(device)
    )
    return model


def embed_sequences(sequences, model_path, device, batch_size=16):
    """
    Embed a list of sequences using a frozen ESM2 model.
    Returns a list of np.arrays of shape (L, hidden_dim), one per sequence.
    Uses SQL mode (same as ultrafast.py) to avoid macOS multiprocessing issues.
    """
    model = load_esm_model(model_path, device)

    db_path = Path(tempfile.NamedTemporaryFile(suffix=".db", delete=False).name)
    model.embed_dataset(
        sequences=sequences,
        tokenizer=model.tokenizer,
        batch_size=batch_size,
        max_len=None,
        full_embeddings=True,
        embed_dtype=torch.float32,
        num_workers=0,
        sql=True,
        sql_db_path=str(db_path),
        save=False,
    )

    conn = sqlite3.connect(db_path)
    embeddings = []
    for seq in sequences:
        row = conn.execute(
            "SELECT embedding, shape, dtype FROM embeddings WHERE sequence = ?", (seq,)
        ).fetchone()
        blob, shape_str, dtype_str = row
        shape = tuple(int(x) for x in shape_str.strip("()").split(",") if x.strip())
        emb = np.frombuffer(blob, dtype=np.dtype(dtype_str)).reshape(shape).copy()
        if emb.shape[0] == len(seq) + 2:
            emb = emb[1:-1]
        embeddings.append(emb)
    conn.close()
    db_path.unlink(missing_ok=True)

    return embeddings


def make_windows(matrix, window_size=11):
    """
    Convert a per-residue embedding matrix into sliding window feature vectors.
    matrix: np.array of shape (L, hidden_dim)
    returns: np.array of shape (L, window_size * hidden_dim)
    """
    pad = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    return np.stack([padded[i : i + window_size].flatten() for i in range(len(matrix))])
