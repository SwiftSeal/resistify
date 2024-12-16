import argparse
from rich_argparse import RichHelpFormatter
import logging
from rich.logging import RichHandler
import sys
import os
from resistify.utility import (
    create_output_directory,
    parse_fasta,
    save_fasta,
    result_table,
    annotation_table,
    domain_table,
    motif_table,
    extract_nbarc,
    coconat_table,
    download_files,
    verify_files,
)
from resistify.hmmsearch import hmmsearch
from resistify.nlrexpress import nlrexpress
from resistify.coconat import coconat
from resistify.tmbed import tmbed
from resistify.__version__ import __version__


def add_common_args(parser):
    """
    Add common arguments shared between NLR and PRR parsers.
    """
    parser.add_argument(
        "input",
        type=validate_input_file,
        help="Path to the input FASTA file containing sequences to be analyzed.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Path to the output directory where results will be saved.",
        default=os.getcwd(),
    )
    parser.add_argument(
        "--models-path",
        help="Path to the downloaded models directory. If not provided, models will be downloaded to $HOME/.cache/ by default.",
        default=None,
    )
    parser.add_argument(
        "--lrr_gap",
        help="Minimum gap (in amino acids) between LRR motifs. Default is 75.",
        default=75,
        type=int,
    )
    parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs required to be considered an LRR domain. Default is 4.",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--duplicate_gap",
        help="Gap size (in amino acids) to consider merging duplicate annotations. Default is 100.",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer.",
        default=None,
        type=int,
    )


def validate_input_file(filepath):
    """
    Validate that the input file exists.
    """
    if not os.path.isfile(filepath):
        raise argparse.ArgumentTypeError(f"Input file {filepath} does not exist")
    return filepath


def parse_args(args=None):
    """
    Parse command-line arguments for Resistify.
    """
    parser = argparse.ArgumentParser(
        description="A tool for identifying and classifying resistance genes in plant genomes.",
        formatter_class=RichHelpFormatter,
    )

    # Global arguments
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Show the version number and exit.",
    )
    parser.add_argument(
        "--debug", help="Enable debug logging for detailed output.", action="store_true"
    )

    # Subparsers
    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Subcommands"
    )

    # NLR subparser
    nlr_parser = subparsers.add_parser(
        "nlr",
        help="Identify and classify NLR resistance genes.",
        formatter_class=RichHelpFormatter,
    )
    nlr_parser.add_argument(
        "--retain",
        help="Non-NLRs will be retained for motif prediction and reported in the final output.",
        action="store_true",
    )
    nlr_parser.add_argument(
        "--coconat",
        help="If enabled, Coconat will be used to improve coiled-coil (CC) annotations.",
        action="store_true",
    )
    nlr_parser.add_argument(
        "--batch",
        help="Number of sequences to process in parallel with coconat.",
        default=None,
        type=int,
    )
    add_common_args(nlr_parser)

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes.",
        formatter_class=RichHelpFormatter,
    )

    # Download models subparser
    download_parser = subparsers.add_parser(
        "download_models",
        help="Download models for CoCoNat and TMbed.",
        formatter_class=RichHelpFormatter,
    )
    download_parser.add_argument(
        "models_path",
        help="Path to the directory which will be used to store downloaded models.",
    )
    add_common_args(prr_parser)

    return parser.parse_args(args)


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


def nlr(args, log):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    sequences = parse_fasta(args.input)
    log.info("Searching for NLRs...")
    sequences = hmmsearch(sequences, "nlr")
    if not args.retain:
        sequences = [sequence for sequence in sequences if sequence.has_nbarc]
        if len(sequences) == 0:
            log.warning("No NLRs detected! Maybe try --retain?")
            return sequences
        else:
            log.info(f"{len(sequences)} NLRs identified...")
    else:
        log.info("NLRexpress will be run against all input sequences...")

    # If not specified, run NLRexpress in batches of 5 on all sequences or 1 on retained NLRs
    if args.retain and args.chunksize is None:
        chunksize = 5
    elif args.chunksize is None:
        chunksize = 1
    else:
        chunksize = args.chunksize

    sequences = nlrexpress(sequences, "all", chunksize)

    if args.coconat:
        log.info("Running CoCoNat to identify additional CC domains...")
        sequences = coconat(sequences, args.models_path)

    log.info("Classifying NLRs...")
    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        if args.coconat:
            sequence.identify_cc_domains()
        sequence.merge_annotations(args.duplicate_gap)
        sequence.classify_nlr()

    return sequences


def prr(args, log):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    if args.chunksize is None:
        chunksize = 5
    else:
        chunksize = args.chunksize

    log.info("Searching for PRRs...")
    sequences = parse_fasta(args.input)
    sequences = hmmsearch(sequences, "prr")

    sequences = tmbed(sequences, args.models_path)
    sequences = [sequence for sequence in sequences if sequence.is_rlp()]
    if len(sequences) > 0:
        log.info(f"{len(sequences)} PRRs identified...")
        sequences = nlrexpress(sequences, "lrr", chunksize)

        log.info("Classifying PRRs...")
        for sequence in sequences:
            sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
            sequence.merge_annotations(args.duplicate_gap)
            sequence.classify_rlp()
    else:
        log.warning("No PRRs detected!")

    return sequences


def download(args, log):
    log.info("Downloading model data...")
    download_files(args.models_path)
    verify_files(args.models_path)
    log.info(
        "Models downloaded successfully. You can supply these to Resistify with the argument `--models <path-to-directory>`"
    )


def main():
    args = parse_args()

    logging.basicConfig(
        level="DEBUG" if args.debug else "INFO",
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(show_path=False)],
    )
    log = logging.getLogger("rich")
    log.info(f"Welcome to Resistify version {__version__}!")

    if args.command == "nlr":
        sequences = nlr(args, log)
    elif args.command == "prr":
        sequences = prr(args, log)
    elif args.command == "download_models":
        download(args, log)
        sys.exit(0)
    else:
        log.error(f"Unknown command: {args.command}")
        sys.exit(1)

    write_results(sequences, args)

    log.info("Thank you for using Resistify!")
    log.info("If you used Resistify in your research, please cite the following:")
    log.info(" - Resistify: https://doi.org/10.1101/2024.02.14.580321")
    log.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if args.command == "prr":
        log.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif args.coconat:
        log.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")


if __name__ == "__main__":
    main()
