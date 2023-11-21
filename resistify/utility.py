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
    parser.add_argument("--evalue", help="E-value threshold for hmmsearch. Scientific notation not accepted!", default="0.00001")
    parser.add_argument("--lrr_gap", help="Gap size for LRR annotation", default=75)
    parser.add_argument("--lrr_length", help="Minimum number of LRR motifs to be considered an LRR domain", default=4)
    parser.add_argument("--duplicate_gap", help="Gap size (aa) to consider merging duplicate annotations", default=100)
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


def save_fasta(sequences, path):
    with open(path, "w") as file:
        for sequence in sequences:
            file.write(f">{sequence}\n")
            file.write(f"{sequences[sequence].sequence}\n")
    return path


def result_table(sequences, results_dir):
    with open(os.path.join(results_dir, "results.tsv"), "w") as file:
        file.write(
            "Sequence\tLength\tDomains\tClassification\tNBARC_motifs\tMADA\tCJID\n"
        )
        for sequence in sequences:
            nbarc_motifs = [
                "VG",
                "P-loop",
                "RNSB-A",
                "Walker-B",
                "RNSB-B",
                "RNSB-C",
                "RNSB-D",
                "GLPL",
                "MHD",
            ]

            n_nbarc_motifs = 0
            for motif in nbarc_motifs:
                if len(sequences[sequence].motifs[motif]) > 0:
                    n_nbarc_motifs += 1

            length = len(sequences[sequence].sequence)
            classification = sequences[sequence].classification
            mada = sequences[sequence].mada
            cjid = sequences[sequence].cjid
            domain_string = sequences[sequence].domain_string

            file.write(
                f"{sequence}\t{length}\t{domain_string}\t{classification}\t{n_nbarc_motifs}\t{mada}\t{cjid}\n"
            )


def domain_table(sequences, results_dir):
    with open(os.path.join(results_dir, "domains.tsv"), "w") as file:
        file.write("Sequence\tDomain\tStart\tEnd\n")
        for sequence in sequences:
            for annotation in sequences[sequence].annotations:
                file.write(
                    f"{sequence}\t{annotation.domain}\t{annotation.start}\t{annotation.end}\n"
                )


def motif_table(sequences, results_dir):
    with open(os.path.join(results_dir, "motifs.tsv"), "w") as file:
        file.write("Sequence\tMotif\tPosition\tProbability\n")
        for sequence in sequences:
            for motif in sequences[sequence].motifs:
                for item in sequences[sequence].motifs[motif]:
                    file.write(
                        f"{sequence}\t{motif}\t{item.position}\t{item.probability}\n"
                    )


def extract_nbarc(sequences, results_dir):
    """
    Extract all nbarc domains of all proteins into a fasta file.
    """
    with open(os.path.join(results_dir, "nbarc.fasta"), "w") as file:
        for sequence in sequences:
            count = 1
            for annotation in sequences[sequence].annotations:
                if annotation.domain == "NB-ARC":
                    file.write(f">{sequence}_{count}\n")
                    file.write(
                        f"{sequences[sequence].sequence[annotation.start-1:annotation.end]}\n"
                    )
                    count += 1
