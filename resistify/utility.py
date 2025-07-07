import sys
from pathlib import Path
import os
import csv
import gzip
from resistify.annotations import Sequence, NBARC_MOTIFS
from resistify._loguru import logger
from resistify.__version__ import __version__


class ProgressLogger:
    def __init__(self, total_count):
        self.total_count = total_count
        self.current_count = 0
        self.last_reported_percent = -1  # Initialize with an invalid percentage

    def update(self):
        self.current_count += 1
        if self.total_count < 10:
            # For small totals, report as "n of total"
            logger.info(f"{self.current_count} of {self.total_count} complete")
        else:
            # Calculate percentage
            percent_complete = int((self.current_count / self.total_count) * 100)
            if self.current_count == self.total_count:
                logger.info("100% complete")
                self.last_reported_percent = 100
            elif (
                percent_complete % 10 == 0
                and percent_complete > 0
                and percent_complete > self.last_reported_percent
            ):
                logger.info(f"{percent_complete}% complete")
                self.last_reported_percent = percent_complete


def hello(models: bool = False):
    """
    Print some introductory information for the user.
    """
    logger.info(f"Welcome to Resistify {__version__}!")
    logger.info("Need help? Visit https://github.com/SwiftSeal/Resistify")

    if models:
        hf_home = os.environ.get("HF_HOME", "~/.cache/huggingface/")
        torch_home = os.environ.get("TORCH_HOME", "~/.cache/torch/")

        logger.info("Models will be downloaded/retrieved from the following locations:")
        logger.info(f"$HF_HOME: {hf_home}")
        logger.info(f"$TORCH_HOME: {torch_home}")


def goodbye(coconat: bool = False, tmbed: bool = False):
    logger.info("Thank you for using Resistify!")
    logger.info("If you used Resistify in your research, please cite the following:")
    logger.info(" - Resistify: https://doi.org/10.1177/11779322241308944")
    logger.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if tmbed:
        logger.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif coconat:
        logger.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")


def create_output_directory(outdir: Path):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        logger.debug(f"Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        logger.error(f"Error creating output directory: {e}")
        sys.exit(1)


def write_results(
    sequences: list,
    outdir: Path,
    command: str,
    coconat: bool = False,
    retain: bool = False,
):
    results_dir = create_output_directory(outdir)
    if command == "nlr":
        result_table(sequences, results_dir, "nlr", retain=retain)
        save_fasta(
            sequences,
            os.path.join(results_dir, "nlr.fasta"),
            classified_only=retain,
        )
        extract_nbarc(sequences, results_dir)
        if coconat:
            coconat_table(sequences, results_dir)
    elif command == "prr":
        result_table(sequences, results_dir, "prr")
        save_fasta(
            sequences, os.path.join(results_dir, "prr.fasta"), classified_only=True
        )
    annotation_table(sequences, results_dir)
    domain_table(sequences, results_dir)
    motif_table(sequences, results_dir)


def check_sequence(sequence):
    if "*" in sequence.seq:
        logger.warning(
            f"An internal '*' character is present in {sequence.id} - skipping this sequence..."
        )
    elif "." in sequence.seq:
        logger.warning(
            f"An internal '.' character is present in {sequence.id} - skipping this sequence..."
        )
    elif len(sequence.seq) > 100000:
        logger.warning(
            f"Sequence {sequence.id} length is > 100,000 - skipping this sequence..."
        )
    elif len(sequence.seq) < 28:
        logger.warning(
            f"Sequence {sequence.id} length is < 28 - skipping this sequence..."
        )
    else:
        return True
    return False


def parse_fasta(infile: Path):
    sequences = []

    # Automatically detect if the file is gzipped
    open_func = gzip.open if infile.suffix == ".gz" else open

    with open_func(infile, "rt") as file:
        seq_id = None
        seq = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    seq = seq.strip("*").strip(".")
                    sequence = Sequence(seq_id, seq)
                    if check_sequence(sequence):
                        sequences.append(Sequence(seq_id, seq))
                seq_id = line[1:].split()[0]
                seq = ""
            else:
                seq += line
        if seq_id:
            seq = seq.strip("*").strip(".")
            sequence = Sequence(seq_id, seq)
            if check_sequence(sequence):
                sequences.append(Sequence(seq_id, seq))

    if not sequences:
        logger.error("No valid sequences found in input file!")
        sys.exit(1)

    return sequences


def wrap_sequence(sequence, wrap_length=80):
    wrapped_sequence = ""
    for i in range(0, len(sequence), wrap_length):
        wrapped_sequence += sequence[i : i + wrap_length] + "\n"
    return wrapped_sequence


def save_fasta(sequences, path, classified_only=False):
    with open(path, "w") as file:
        for sequence in sequences:
            # Special case for PRRs as we are interested in type
            if classified_only and sequence.type in ["RLP", "RLK"]:
                file.write(f">{sequence.id}\n")
                file.write(f"{wrap_sequence(sequence.seq)}")
            elif classified_only and sequence.classification is None:
                continue
            else:
                file.write(f">{sequence.id}\n")
                file.write(f"{wrap_sequence(sequence.seq)}")
    return path


def result_table(sequences, results_dir, type, retain=False):
    with open(os.path.join(results_dir, "results.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        if type == "nlr":
            table_writer.writerow(
                [
                    "Sequence",
                    "Length",
                    "LRR_Length",
                    "Motifs",
                    "Domains",
                    "Classification",
                    "NBARC_motifs",
                    "MADA",
                    "MADAL",
                    "CJID",
                ]
            )

            for sequence in sequences:
                if sequence.type == "NLR" or retain:
                    n_nbarc_motifs = 0
                    for motif in NBARC_MOTIFS:
                        if len(sequence.motifs[motif]) > 0:
                            n_nbarc_motifs += 1

                    table_writer.writerow(
                        [
                            sequence.id,
                            len(sequence.seq),
                            sequence.lrr_length,
                            sequence.motif_string,
                            sequence.domain_string,
                            sequence.classification,
                            n_nbarc_motifs,
                            sequence.has_mada,
                            sequence.has_madal,
                            sequence.has_cjid,
                        ]
                    )
        elif type == "prr":
            table_writer.writerow(
                [
                    "Sequence",
                    "Length",
                    "Extracellular_Length",
                    "LRR_Length",
                    "Type",
                    "Classification",
                    "Signal_peptide",
                ]
            )
            for sequence in sequences:
                if sequence.type == "RLP" or sequence.type == "RLK":
                    table_writer.writerow(
                        [
                            sequence.id,
                            len(sequence.seq),
                            sequence.extracellular_length,
                            sequence.lrr_length,
                            sequence.type,
                            sequence.classification,
                            sequence.has_signal_peptide,
                        ]
                    )


def domain_table(sequences, results_dir):
    with open(os.path.join(results_dir, "domains.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(["Sequence", "Domain", "Start", "End"])
        for sequence in sequences:
            for annotation in sequence.merged_annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                    ]
                )


def annotation_table(sequences, results_dir):
    with open(os.path.join(results_dir, "annotations.tsv"), "w") as file:
        table_writer = csv.writer(file, delimiter="\t")
        table_writer.writerow(
            ["Sequence", "Domain", "Start", "End", "E_value", "Score", "Source"]
        )
        for sequence in sequences:
            for annotation in sequence.annotations:
                table_writer.writerow(
                    [
                        sequence.id,
                        annotation.domain,
                        annotation.start,
                        annotation.end,
                        annotation.evalue,
                        annotation.score,
                        annotation.source,
                    ]
                )


def motif_table(sequences, results_dir):
    from resistify.nlrexpress import MOTIF_SPAN_LENGTHS

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
            for motif in sequence.motifs:
                for item in sequence.motifs[motif]:
                    aa_sequence = sequence.seq
                    downstream_sequence = aa_sequence[item.position - 5 : item.position]
                    motif_sequence = aa_sequence[
                        item.position : item.position + MOTIF_SPAN_LENGTHS[motif]
                    ]
                    upstream_sequence = aa_sequence[
                        item.position + MOTIF_SPAN_LENGTHS[motif] : item.position
                        + MOTIF_SPAN_LENGTHS[motif]
                        + 5
                    ]
                    table_writer.writerow(
                        [
                            sequence.id,
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
            for annotation in sequence.merged_annotations:
                if annotation.domain == "NB-ARC":
                    file.write(f">{sequence.id}_{count}\n")
                    file.write(
                        f"{wrap_sequence(sequence.seq[annotation.start : annotation.end])}"
                    )
                    count += 1


def coconat_table(sequences, results_dir):
    output_path = os.path.join(results_dir, "coconat.tsv")
    with open(output_path, "w") as f:
        f.write("Sequence\tPosition\tProbability\n")
        for sequence in sequences:
            for i, probability in enumerate(sequence.cc_probs):
                f.write(f"{sequence.id}\t{i}\t{probability}\n")
