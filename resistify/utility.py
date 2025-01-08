import sys
import os
import csv
import json
import hashlib
import requests
import gzip
from resistify.annotations import Sequence, NBARC_MOTIFS
from resistify._loguru import logger


def logger_format(debug):
    logging.basicConfig(
        level="DEBUG" if debug else "INFO",
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="[%H:%M:%S]",
    )


def log_percentage(n, total):
    if total < 10:
        logger.info(f"{n} of {total} complete")
    elif n % (total // 10) == 0:
        percent_complete = n / total * 100
        logger.info(f"{int(round(percent_complete, -1))}% complete")


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        logger.debug(f"Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        logger.error(f"Error creating output directory: {e}")
        sys.exit(1)


def write_results(sequences, args):
    results_dir = create_output_directory(args.outdir)
    if args.command == "nlr":
        result_table(sequences, results_dir, "nlr", retain=args.retain)
        save_fasta(
            sequences,
            os.path.join(results_dir, "nlr.fasta"),
            classified_only=args.retain,
        )
        extract_nbarc(sequences, results_dir)
        if args.coconat:
            coconat_table(sequences, results_dir)
    elif args.command == "prr":
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


def parse_fasta(path):
    sequences = []

    # Automatically detect if the file is gzipped
    open_func = gzip.open if path.endswith(".gz") else open

    with open_func(path, "rt") as file:
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
                        f"{wrap_sequence(sequence.seq[annotation.start:annotation.end])}"
                    )
                    count += 1


def coconat_table(sequences, results_dir):
    output_path = os.path.join(results_dir, "coconat.tsv")
    with open(output_path, "w") as f:
        f.write("Sequence\tPosition\tProbability\n")
        for sequence in sequences:
            for i, probability in enumerate(sequence.cc_probs):
                f.write(f"{sequence.id}\t{i}\t{probability}\n")


def download_files(base_download_path):
    config_path = os.path.join(os.path.dirname(__file__), "data", "model_paths.json")

    # Read the configuration
    with open(config_path, "r") as f:
        config = json.load(f)

    # Ensure base download path exists
    os.makedirs(base_download_path, exist_ok=True)

    # Process each model's files
    for model_name, model_config in config.items():
        # Create model-specific directory
        model_dir = os.path.join(base_download_path, model_config["directory"])
        os.makedirs(model_dir, exist_ok=True)

        # Download each file
        for file_info in model_config["files"]:
            url = file_info["url"]
            filename = os.path.basename(url)
            local_path = os.path.join(model_dir, filename)

            # Download file
            response = requests.get(url, stream=True)
            response.raise_for_status()

            with open(local_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)

            # Calculate SHA256
            with open(local_path, "rb") as f:
                file_hash = hashlib.sha256()
                for chunk in iter(lambda: f.read(4096), b""):
                    file_hash.update(chunk)

            # Update file info with local path and hash
            file_info["local_path"] = local_path
            file_info["sha256"] = file_hash.hexdigest()

    return config


def verify_files(base_download_path: str):
    # Initialize results
    config_path = os.path.join(os.path.dirname(__file__), "data", "model_paths.json")

    with open(config_path, "r") as f:
        config = json.load(f)

    # Check each model's files
    for model_name, model_config in config.items():
        model_dir = os.path.join(base_download_path, model_config["directory"])

        # Check if model directory exists
        if not os.path.exists(model_dir):
            logger.error(f"Directory not found for {model_name}: {model_dir}")
            sys.exit(1)

        # Check each file
        for file_info in model_config["files"]:
            url = file_info["url"]
            filename = os.path.basename(url)
            local_path = os.path.join(model_dir, filename)

            # Check if file exists
            if not os.path.exists(local_path):
                logger.error(f"File not found: {local_path}")
                sys.exit(1)

    return
