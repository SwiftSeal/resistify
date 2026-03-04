"""
Convert NLRexpress .input files to per-motif .labels files.
Reads training_sets/{cc,lrr,nbs}_e5.input and writes training/data/<motif>.labels
"""

from collections import defaultdict
from pathlib import Path

OUTPUT_DIR = Path("training/data")

# Each input file has its own label map — files must be kept separate since
# sequences from one domain class should not appear as negatives in another.
FILE_CONFIGS = {
    Path("training_sets/cc_e5.input"): {
        "rdhhhdhEDVID": "extEDVID",
    },
    Path("training_sets/nbs_e5.input"): {
        "GmGGvGKTT": "P-loop",
        "GLPLA": "GLPL",
        "MHD": "MHD",
        "KRhhhhDD": "Walker-B",
        "hhGRE": "VG",
        "FDhrhWhshs": "RNBS-A",
        "KhhhTTR": "RNBS-B",
        "LseeeSWeLF": "RNBS-C",
        "CFLYCSLFP": "RNBS-D",
    },
    Path("training_sets/lrr_e5.input"): {
        "L": "LxxLxL",
        "N": "LxxLxL",
        "C": "LxxLxL",
    },
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for input_file, label_map in FILE_CONFIGS.items():
    print(f"Reading {input_file}...")
    motif_names = set(label_map.values())

    # labels[motif][seq_id] = [(pos, residue, label), ...]
    labels = {motif: defaultdict(list) for motif in motif_names}

    with open(input_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            seq_id, pos, residue, raw_label = parts[0], parts[1], parts[2], parts[3]
            for motif in motif_names:
                is_positive = raw_label in label_map and label_map[raw_label] == motif
                labels[motif][seq_id].append((pos, residue, 1 if is_positive else 0))

    for motif, motif_labels in labels.items():
        labels_path = OUTPUT_DIR / f"{motif}.labels"
        n_pos = sum(r[2] for rows in motif_labels.values() for r in rows)
        with open(labels_path, "w") as f:
            f.write("seq_id\tpos\tresidue\tlabel\n")
            for seq_id, rows in motif_labels.items():
                for pos, residue, label in rows:
                    f.write(f"{seq_id}\t{pos}\t{residue}\t{label}\n")
        print(f"  {motif}: {n_pos} positive -> {labels_path}")
