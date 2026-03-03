"""
Data loading utilities for motif classifier training.
"""
from collections import defaultdict


def parse_fasta(path):
    """Parse a FASTA file into a dict of {seq_id: sequence}."""
    seqs, current_id = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:].split()[0]
                seqs[current_id] = []
            elif current_id:
                seqs[current_id].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


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
            seq_id, pos, residue, label = parts[0], int(parts[1]), parts[2], parts[3]
            data[seq_id].append((pos, residue, int(parts[3])))
    return dict(data)