import os
import platform
import torch
import argparse
from pathlib import Path
import logging
from resistify.__version__ import __version__
from resistify.annotation import save_results
from resistify.hmmer import hmmsearch
from resistify.nlrexpress import nlrexpress
from resistify.coconat import predict_coils
from resistify.tmbed import tmbed
from resistify.hmmer import NLR_HMM_DB, RLP_HMM_DB

# rich is optional
try:
    from resistify.console import console
    from rich.logging import RichHandler
    from rich_argparse import ArgumentDefaultsRichHelpFormatter

    help_formatter = ArgumentDefaultsRichHelpFormatter
    logging_handler = [RichHandler(console=console)]
except ImportError:
    help_formatter = argparse.ArgumentDefaultsHelpFormatter
    logging_handler = None

logger = logging.getLogger(__name__)

DEFAULT_DEVICE = str(
    torch.device(
        "mps"
        if torch.backends.mps.is_built()
        else "cuda"
        if torch.cuda.is_available()
        else "cpu"
    )
)


def add_common_args(parser):
    """
    Add common arguments shared between NLR and PRR parsers.
    """
    parser.add_argument(
        "input",
        type=Path,
        help="Path to the input FASTA file containing sequences to be analysed",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        help="Path to the output directory where results will be saved",
        default=Path(os.getcwd()),
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="Number of threads to use for HMMER searches. Default is 0 (auto-detect)",
        default=0,
    )
    parser.add_argument(
        "--device",
        type=str,
        default=DEFAULT_DEVICE,
        help="Torch device to be used - can be 'cpu', 'cuda', or 'mps'",
    )
    parser.add_argument(
        "--lrr_gap",
        help="Minimum gap (in amino acids) between LRR motifs. Default is 75",
        default=75,
        type=int,
    )
    parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs required to be considered an LRR domain. Default is 4",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--duplicate_gap",
        help="Gap size (in amino acids) to consider merging duplicate annotations. Default is 100",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--pfam",
        type=Path,
        default=None,
        help="Path to the Pfam database file - by default a smaller task-specific database is used",
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="A tool for identifying and classifying resistance genes in plant genomes",
        formatter_class=help_formatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"v{__version__}",
        help="Show the version number and exit",
    )
    parser.add_argument(
        "--debug", help="Enable debug logging for detailed output.", action="store_true"
    )

    subparsers = parser.add_subparsers(
        dest="command", required=True, help="Subcommands"
    )

    # NLR subparser
    nlr_parser = subparsers.add_parser(
        "nlr",
        help="Identify and classify NLR resistance genes",
        formatter_class=parser.formatter_class,
    )
    nlr_parser.add_argument(
        "--retain",
        help="Non-NLRs will be retained for motif prediction and reported in the final output.",
        action="store_true",
    )
    nlr_parser.add_argument(
        "--coconat",
        action="store_true",
        help="Enable CoCoNat coiled-coil (CC) domain predictions. CCs will be predicted in the N-terminal region only",
    )
    add_common_args(nlr_parser)

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes",
        formatter_class=parser.formatter_class,
    )
    add_common_args(prr_parser)

    return parser.parse_args()


def main():
    args = parse_args()

    level = logging.DEBUG if args.debug else logging.INFO
    FORMAT = "%(message)s"
    logging.basicConfig(
        level=level, format=FORMAT, datefmt="[%X]", handlers=logging_handler
    )

    logger.info("Resistify 2.0.0")
    logger.info(f"Python version: {platform.python_version()}")
    logger.info(f"PyTorch version: {torch.__version__}")
    logger.info(f"Using device: {args.device}")
    logger.info(f"OS: {platform.system()} {platform.machine()}")

    if args.command == "nlr":
        hmmer_db = args.pfam if args.pfam is not None else NLR_HMM_DB
        proteins = hmmsearch(args.input, hmmer_db, threads=args.threads)

        if not args.retain:
            proteins = {k: v for k, v in proteins.items() if v.is_nlr()}

        nlrexpress(proteins, threads=args.threads)

        if args.coconat:
            predict_coils(proteins, args.device)

        for protein in proteins.values():
            protein.annotate_lrr(lrr_gap=args.lrr_gap, lrr_length=args.lrr_length)
            protein.merge_domains()
            protein.classify_nlr()

    elif args.command == "prr":
        hmmer_db = args.pfam if args.pfam is not None else RLP_HMM_DB
        proteins = hmmsearch(args.input, hmmer_db, threads=args.threads)
        tmbed(proteins, device=args.device)
        proteins = {k: v for k, v in proteins.items() if v.is_rlp()}

        nlrexpress(proteins, search_type="lrr", threads=args.threads)
        for protein in proteins.values():
            protein.annotate_lrr(lrr_gap=args.lrr_gap, lrr_length=args.lrr_length)
            protein.merge_domains()
            protein.classify_rlp()
    else:
        return

    save_results(proteins, args.outdir, command=args.command)

    logger.info("Thank you for using Resistify!")
    logger.info("If you used Resistify in your research, please cite the following:")
    logger.info(" - Resistify: https://doi.org/10.1177/11779322241308944")
    logger.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if args.command == "prr":
        logger.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif args.coconat:
        logger.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")
