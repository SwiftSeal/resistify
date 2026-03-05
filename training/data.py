"""
Data loading utilities for motif classifier training.
"""

from collections import defaultdict


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
