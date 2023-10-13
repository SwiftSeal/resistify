import logging
import argparse
import sys
import os
from Bio import SeqIO
from resistify.annotations import Sequence

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for predicting resistance genes in plant genomes.
        Domains are predicted using hmmsearch against multiple databases which are classified as common resistance gene related domains.
        Sequences are then classified based on the presence and order of these domains.
        """
    )
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("outdir", help="Output directory")

    return parser.parse_args()


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        logging.info(f"😊 Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        logging.error(f"😞 Error creating output directory: {e}")
        sys.exit(1)


def parse_fasta(path):
    sequences = {}
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences[record.id] = Sequence(record.seq)
    return sequences


def save_fasta(sequences, path, subset=None):
    """
    Save a dictionary of sequences to a FASTA file.
    Can optionally save only a subset of the sequences.
    """
    if subset is None:
        with open(path, "w") as file:
            for sequence in sequences:
                file.write(f">{sequence}\n")
                file.write(f"{sequences[sequence].sequence}\n")
    else:
        subset_sequences = {key: sequences[key] for key in subset}
        with open(path, "w") as file:
            for sequence in subset_sequences:
                file.write(f">{sequence}\n")
                file.write(f"{subset_sequences[sequence].sequence}\n")
    return path

