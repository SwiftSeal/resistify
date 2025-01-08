import argparse
import sys
import os
from resistify.utility import (
    write_results,
    parse_fasta,
    download_files,
    verify_files,
)
from resistify._loguru import logger
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
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads available for nlrexpress. Default is the number of available CPUs.",
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
    add_common_args(nlr_parser)

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes.",
    )

    # Download models subparser
    download_parser = subparsers.add_parser(
        "download_models",
        help="Download models for CoCoNat and TMbed.",
    )
    download_parser.add_argument(
        "models_path",
        help="Path to the directory which will be used to store downloaded models.",
    )
    add_common_args(prr_parser)

    return parser.parse_args(args)

def nlr(args):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    sequences = parse_fasta(args.input)
    logger.info("Searching for NLRs...")
    sequences = hmmsearch(sequences, "nlr")
    if not args.retain:
        sequences = [sequence for sequence in sequences if sequence.has_nbarc]
        if len(sequences) == 0:
            logger.warning("No NLRs detected! Maybe try --retain?")
            return sequences
        else:
            logger.info(f"{len(sequences)} NLRs identified...")
    else:
        logger.info("NLRexpress will be run against all input sequences...")

    # If not specified, run NLRexpress in batches of 5 on all sequences or 1 on retained NLRs
    if args.retain and args.chunksize is None:
        chunksize = 5
    elif args.chunksize is None:
        chunksize = 1
    else:
        chunksize = args.chunksize

    sequences = nlrexpress(sequences, "all", chunksize, args.threads, args.debug)

    if args.coconat:
        logger.info("Running CoCoNat to identify additional CC domains...")
        sequences = coconat(sequences, args.models_path)

    logger.info("Classifying NLRs...")
    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        if args.coconat:
            sequence.identify_cc_domains()
        sequence.merge_annotations(args.duplicate_gap)
        sequence.classify_nlr()

    return sequences


def prr(args):
    # Check to see if all provided models exist
    if args.models_path is not None:
        verify_files(args.models_path)

    if args.chunksize is None:
        chunksize = 5
    else:
        chunksize = args.chunksize

    logger.info("Searching for PRRs...")
    sequences = parse_fasta(args.input)
    sequences = hmmsearch(sequences, "prr")

    sequences = tmbed(sequences, args.models_path)
    sequences = [sequence for sequence in sequences if sequence.is_rlp()]
    if len(sequences) > 0:
        logger.info(f"{len(sequences)} PRRs identified...")
        sequences = nlrexpress(sequences, "lrr", chunksize, args.threads, args.debug)

        logger.info("Classifying PRRs...")
        for sequence in sequences:
            sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
            sequence.merge_annotations(args.duplicate_gap)
            sequence.classify_rlp()
    else:
        logger.warning("No PRRs detected!")

    return sequences


def download(args):
    logger.info("Downloading model data...")
    download_files(args.models_path)
    verify_files(args.models_path)
    logger.info(
        "Models downloaded successfully. You can supply these to Resistify with the argument `--models <path-to-directory>`"
    )


def main():
    args = parse_args()

    if args.debug:
        logger.update_level("DEBUG")

    logger.info(f"Welcome to Resistify version {__version__}!")

    if args.command == "nlr":
        sequences = nlr(args)
    elif args.command == "prr":
        sequences = prr(args)
    elif args.command == "download_models":
        download(args)
        sys.exit(0)
    else:
        logger.error(f"Unknown command: {args.command}")
        sys.exit(1)

    write_results(sequences, args)

    logger.info("Thank you for using Resistify!")
    logger.info("If you used Resistify in your research, please cite the following:")
    logger.info(" - Resistify: https://doi.org/10.1101/2024.02.14.580321")
    logger.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if args.command == "prr":
        logger.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif args.coconat:
        logger.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")


if __name__ == "__main__":
    main()
