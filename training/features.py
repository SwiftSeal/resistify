"""
Feature extraction utilities for motif classifier training.
"""
import numpy as np
import torch
from transformers import AutoModel


def load_esm_model(model_path, device):
    """Load a frozen ESM2 model onto the specified device."""
    print(f"Loading {model_path} on {device}...")
    model = AutoModel.from_pretrained(
        model_path,
        trust_remote_code=True,
    ).half().eval().to(device)
    return model


def embed_sequences(sequences, model_path, device, batch_size=16):
    """
    Embed a list of sequences using a frozen ESM2 model.
    Returns a list of np.arrays of shape (L, hidden_dim), one per sequence.
    """
    model     = load_esm_model(model_path, device)
    tokenizer = model.tokenizer

    result = model.embed_dataset(
        sequences       = sequences,
        tokenizer       = tokenizer,
        batch_size      = batch_size,
        max_len         = None,
        full_embeddings = True,
        embed_dtype     = torch.float32,
        num_workers     = 0,
        save            = False,
    )

    # embed_dataset sorts sequences internally - restore original order
    embedding_map = {}
    for seq, emb in result.items():
        emb_np = emb.numpy() if isinstance(emb, torch.Tensor) else np.array(emb)
        if emb_np.shape[0] == len(seq) + 2:
            emb_np = emb_np[1:-1]  # strip BOS/EOS
        embedding_map[seq] = emb_np

    return [embedding_map[s] for s in sequences]


def make_windows(matrix, window_size=11):
    """
    Convert a per-residue embedding matrix into sliding window feature vectors.
    matrix: np.array of shape (L, hidden_dim)
    returns: np.array of shape (L, window_size * hidden_dim)
    """
    pad    = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    return np.stack([
        padded[i:i + window_size].flatten() for i in range(len(matrix))
    ])