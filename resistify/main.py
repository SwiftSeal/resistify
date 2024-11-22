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

__version__ = "0.6.0"


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
        "--evalue",
        help="E-value threshold for hmmsearch. Default is 0.00001.",
        default="0.00001",
    )
    parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer. Default is 5.",
        default=5,
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
    add_common_args(nlr_parser)
    nlr_parser.add_argument(
        "--retain",
        help="Non-NLRs will be retained for motif prediction and reported in the final output.",
        action="store_true",
    )
    nlr_parser.add_argument(
        "--batch",
        help="Number of sequences to process in parallel with jackhmmer.",
        default=None,
        type=int,
    )
    nlr_parser.add_argument(
        "--coconat",
        help="If enabled, Coconat will be used to improve coiled-coil (CC) annotations.",
        action="store_true",
    )

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes.",
        formatter_class=RichHelpFormatter,
    )
    add_common_args(prr_parser)

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

    return parser.parse_args(args)


def nlr(args, log):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    sequences = parse_fasta(args.input)
    log.info("Searching for NLRs...")
    sequences = hmmsearch(sequences, "nlr", args.evalue)
    if not args.retain:
        sequences = [sequence for sequence in sequences if sequence.has_nbarc()]
        if len(sequences) == 0:
            log.error("No NLRs detected! Maybe try --retain?")
            sys.exit(1)
        else:
            log.info(f"{len(sequences)} NLRs identified...")
    else:
        log.info("NLRexpress will be run against all input sequences...")

    sequences = nlrexpress(
        sequences,
        "all",
        args.chunksize,
    )

    if args.coconat:
        log.info("Running CoCoNat to identify additional CC domains...")
        sequences = coconat(sequences, args.models_path)

    log.info("Classifying sequences...")
    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        if args.coconat:
            sequence.identify_cc_domains()
        sequence.merge_annotations(args.duplicate_gap)
        sequence.classify_nlr()

    results_dir = create_output_directory(args.outdir)
    log.info(f"Saving results to {results_dir}")
    result_table(sequences, results_dir, "nlr", retain=args.retain)
    annotation_table(sequences, results_dir)
    domain_table(sequences, results_dir)
    motif_table(sequences, results_dir)
    extract_nbarc(sequences, results_dir)
    save_fasta(sequences, os.path.join(results_dir, "nlr.fasta"), classified_only=True)
    if args.coconat:
        coconat_table(sequences, results_dir)


def prr(args, log):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    log.info("Searching for PRRs...")
    sequences = parse_fasta(args.input)
    sequences = hmmsearch(sequences, "prr", args.evalue)
    sequences = nlrexpress(
        sequences,
        "lrr",
        args.chunksize,
    )

    sequences = tmbed(
        sequences, args.models_path
    )  # Right for some reason if this precedes nlrexpress(), it freezes? dunno why but just make sure it's downstream...

    sequences = [sequence for sequence in sequences if sequence.is_rlp()]

    log.info("Classifying sequences...")
    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        sequence.merge_annotations(args.duplicate_gap)
        sequence.classify_rlp()

    results_dir = create_output_directory(args.outdir)
    log.info(f"Saving results to {results_dir}")
    result_table(sequences, results_dir, "prr")
    annotation_table(sequences, results_dir)
    domain_table(sequences, results_dir)
    motif_table(sequences, results_dir)
    save_fasta(sequences, os.path.join(results_dir, "prr.fasta"), classified_only=True)


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
        handlers=[RichHandler()],
    )
    log = logging.getLogger("rich")
    log.info(f"Welcome to Resistify version {__version__}!")

    if args.command == "nlr":
        nlr(args, log)
    elif args.command == "prr":
        prr(args, log)
    elif args.command == "download_models":
        download(args, log)
        sys.exit(0)

    else:
        log.error(f"Unknown command: {args.command}")
        sys.exit(1)

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
