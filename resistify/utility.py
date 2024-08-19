import sys
import os
import csv
import shutil
import logging
from Bio import SeqIO
from resistify.annotations import Sequence
from resistify.nlrexpress import MOTIF_SPAN_LENGTHS
from tempfile import TemporaryDirectory

log = logging.getLogger(__name__)


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        log.debug(f"Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        log.error(f"Error creating output directory: {e}")
        sys.exit(1)


def prepare_temp_directory(data_dir):
    temp_dir = TemporaryDirectory()

    # Copy nlrexpress database to the temporary directory
    log.debug(f"Copying nlrexpress database to {temp_dir.name}")
    shutil.copy(
        os.path.join(data_dir, "nlrexpress.fasta"),
        os.path.join(temp_dir.name, "nlrexpress.fasta"),
    )
    return temp_dir


def parse_fasta(path):
    sequences = {}
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            # need to remove asterisk, interferes with hmmsearch
            sequence_str = str(record.seq).strip("*")
            if "*" in sequence_str:
                log.error(f"Internal stop codon detected in sequence {record.id}")
                sys.exit(1)
            sequences[record.id] = Sequence(sequence_str)

    if len(sequences) > 100000:
        log.error("Input fasta larger than 100,000 sequences currently not supported.")
        sys.exit(1)

    return sequences


def save_fasta(sequences, path):
    with open(path, "w") as file:
        for sequence in sequences:
            file.write(f">{sequence}\n")
            file.write(f"{sequences[sequence].sequence}\n")
    return path


def result_table(sequences, results_dir):
    with open(os.path.join(results_dir, "results.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            [
                "Sequence",
                "Length",
                "Motifs",
                "Domains",
                "Classification",
                "NBARC_motifs",
                "MADA",
                "MADAL",
                "CJID",
            ]
        )

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

        for sequence in sequences:
            n_nbarc_motifs = 0
            for motif in nbarc_motifs:
                if len(sequences[sequence].motifs[motif]) > 0:
                    n_nbarc_motifs += 1

            length = len(sequences[sequence].sequence)
            motif_string = sequences[sequence].motif_string
            classification = sequences[sequence].classification
            mada = sequences[sequence].mada
            madal = sequences[sequence].madal
            cjid = sequences[sequence].cjid
            domain_string = sequences[sequence].domain_string

            table_writer.writerow(
                [
                    sequence,
                    length,
                    motif_string,
                    domain_string,
                    classification,
                    n_nbarc_motifs,
                    mada,
                    madal,
                    cjid,
                ]
            )


def domain_table(sequences, results_dir):
    with open(os.path.join(results_dir, "domains.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(["Sequence", "Domain", "Start", "End", "E_value"])
        for sequence in sequences:
            for annotation in sequences[sequence].annotations:
                table_writer.writerow(
                    [
                        sequence,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                        annotation.evalue,
                    ]
                )


def motif_table(sequences, results_dir):
    output_path = os.path.join(results_dir, "motifs.tsv")
    with open(output_path, "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            [
                "Sequence",
                "Motif",
                "Position",
                "Probability",
                "Downstream_sequence",
                "Motif_sequence",
                "Upstream_sequence",
            ]
        )
        for sequence in sequences:
            for motif in sequences[sequence].motifs:
                for item in sequences[sequence].motifs[motif]:
                    aa_sequence = sequences[sequence].sequence
                    downstream_sequence = aa_sequence[item.position - 5 : item.position]
                    motif_sequence = aa_sequence[
                        item.position : item.position + MOTIF_SPAN_LENGTHS[motif]
                    ]
                    upstream_sequence = aa_sequence[
                        item.position
                        + MOTIF_SPAN_LENGTHS[motif] : item.position
                        + MOTIF_SPAN_LENGTHS[motif]
                        + 5
                    ]
                    table_writer.writerow(
                        [
                            sequence,
                            motif,
                            item.position,
                            item.probability,
                            downstream_sequence,
                            motif_sequence,
                            upstream_sequence,
                        ]
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
                        f"{sequences[sequence].sequence[annotation.start:annotation.end]}\n"
                    )
                    count += 1
