"""
Data loading utilities for motif classifier training.
"""

from collections import defaultdict
import numpy as np


def parse_labels(path):
    """
    Parse a .labels file into a dict of {seq_id: [(pos, residue, label), ...]}.
    Expects a tab-separated file with header: seq_id, pos, residue, label
    """
    data = defaultdict(list)
    with open(path) as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            seq_id, pos, residue, label = (
                parts[0],
                int(parts[1]),
                parts[2],
                int(parts[3]),
            )
            data[seq_id].append((pos, residue, label))
    return dict(data)


def augment_label_data(
    label_data, rng=None, split_prob=0.5, n_concat_ratio=1.0, min_chunk=20
):
    """
    Augment label_data with two strategies:

    1. Split: each sequence is split into 2–3 random chunks with probability
       `split_prob`.  Each chunk becomes an independent synthetic sequence.
    2. Concat: random pairs of sequences are concatenated to form new sequences.
       The number of concatenations is len(label_data) * n_concat_ratio.

    Returns
    -------
    aug_label_data : dict
        New entries only (originals not included).
    group_map : dict
        Maps every seq_id (original + augmented) to an original seq_id so that
        GroupShuffleSplit keeps related sequences in the same fold.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    aug_data = {}
    # originals map to themselves
    group_map = {s: s for s in label_data}

    seq_ids = list(label_data.keys())

    # ── Split augmentation ────────────────────────────────────────────────────
    for seq_id in seq_ids:
        labels = label_data[seq_id]
        L = len(labels)
        if L < min_chunk * 2:
            continue
        if rng.random() > split_prob:
            continue

        n_chunks = int(rng.choice([2, 3]))
        # ensure we can carve out chunks of at least min_chunk residues each
        if L < min_chunk * n_chunks:
            n_chunks = 2

        if n_chunks == 2:
            sp = int(rng.integers(min_chunk, L - min_chunk))
            splits = [labels[:sp], labels[sp:]]
        else:
            sp1 = int(rng.integers(min_chunk, L - 2 * min_chunk))
            sp2 = int(rng.integers(sp1 + min_chunk, L - min_chunk))
            splits = [labels[:sp1], labels[sp1:sp2], labels[sp2:]]

        for i, chunk in enumerate(splits):
            aug_id = f"{seq_id}__split{i}"
            # re-index positions from 0 so downstream code stays consistent
            aug_data[aug_id] = [(j, r, lbl) for j, (_, r, lbl) in enumerate(chunk)]
            group_map[aug_id] = seq_id

    # ── Concatenation augmentation ────────────────────────────────────────────
    n_concat = max(1, int(len(seq_ids) * n_concat_ratio))
    for _ in range(n_concat):
        id1, id2 = rng.choice(seq_ids, size=2, replace=False)
        labels1 = label_data[id1]
        labels2 = label_data[id2]
        offset = len(labels1)
        combined = labels1 + [(offset + p, r, lbl) for p, r, lbl in labels2]
        aug_id = f"{id1}__cat__{id2}"
        aug_data[aug_id] = combined
        # group by the first (primary) sequence to avoid val leakage from id1
        group_map[aug_id] = id1

    return aug_data, group_map
