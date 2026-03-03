"""
Convert .input files to per-motif .fasta and .labels files.
Usage: python convert_input.py nbs_e5.input training/data/nbs
"""

import sys
from collections import defaultdict

input_file = sys.argv[1]
out_prefix = sys.argv[2]

LABEL_MAP = {
    "rdhhhdhEDVID": "extEDVID",
    "GmGGvGKTT": "P-loop",
    "GLPLA": "GLPL",
    "MHD": "MHD",
    "KRhhhhDD": "Walker-B",
    "hhGRE": "RNSB-A",
    "FDhrhWhshs": "RNSB-B",
    "LseeeSWeLF": "RNSB-C",
    "KhhhTTR": "RNSB-D",
    "CFLYCSLFP": "VG",
    "LxxLxL": "lrr",
}

sequences = defaultdict(list)
labels = defaultdict(lambda: defaultdict(list))

with open(input_file) as f:
    for line in f:
        parts = line.strip().split("\t")
        seq_id, pos, residue, raw_label = parts[0], parts[1], parts[2], parts[3]
        sequences[seq_id].append(residue)

        for motif in LABEL_MAP.values():
            is_positive = raw_label in LABEL_MAP and LABEL_MAP[raw_label] == motif
            labels[motif][seq_id].append((pos, residue, 1 if is_positive else 0))

# Shared FASTA
fasta_path = f"{out_prefix}.fasta"
with open(fasta_path, "w") as f:
    for seq_id, residues in sequences.items():
        f.write(f">{seq_id}\n{''.join(residues)}\n")
print(f"Written {len(sequences)} sequences to {fasta_path}")

# Per-motif label files
for motif, motif_labels in labels.items():
    labels_path = f"{out_prefix}_{motif}.labels"
    n_pos = sum(r[2] for rows in motif_labels.values() for r in rows)
    with open(labels_path, "w") as f:
        f.write("seq_id\tpos\tresidue\tlabel\n")
        for seq_id, rows in motif_labels.items():
            for pos, residue, label in rows:
                f.write(f"{seq_id}\t{pos}\t{residue}\t{label}\n")
    print(f"  {motif}: {n_pos} positive -> {labels_path}")
