import os
import platform
import argparse
from pathlib import Path
import logging
from resistify.__version__ import __version__
from resistify.parse_fasta import parse_fasta
from resistify.annotation import save_results
from resistify.hmmer import hmmsearch
from resistify.nlrexpress import nlrexpress
from resistify.coconat import predict_coils
from resistify.tmbed import tmbed
from resistify.hmmer import NLR_HMM_DB, RLP_HMM_DB
from resistify.device import get_device, get_threads

logger = logging.getLogger(__name__)

def add_common_args(parser):
    """
    Add common arguments shared between NLR and PRR parsers.
    """
    parser.add_argument(
        "input",
        help="Path to the input FASTA file containing sequences to be analysed",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--outdir",
        help="Path to the output directory where results will be saved",
        type=Path,
        default=Path(os.getcwd()),
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Number of threads to use for HMMER searches.",
        type=int,
        default=get_threads(),
    )
    parser.add_argument(
        "--device",
        help="Device to use for CoCoNat and TMbed predictions. Selects the best available device by default.",
        type=str,
        default=get_device(),
        choices=["cpu", "cuda", "mps"],
    )
    parser.add_argument(
        "--batch_size",
        help="Batch size for CoCoNat and TMbed predictions. Adjust based on available GPU memory.",
        type=int,
        default=4,
    )
    parser.add_argument(
        "--lrr_gap",
        help="Minimum gap (in amino acids) between LRR motifs.",
        default=75,
        type=int,
    )
    parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs required to be considered an LRR domain.",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--duplicate_gap",
        help="Gap size (in amino acids) to consider merging duplicate annotations.",
        default=100,
        type=int,
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="A tool for identifying and classifying resistance genes in plant genomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
    FORMAT = "%(asctime)s %(levelname)s - %(message)s"
    logging.basicConfig(level=level, format=FORMAT, datefmt="[%X]")

    logger.info("Welcome to Resistify 2.0.0!")
    logger.info(f"Python version: {platform.python_version()}")
    logger.info(f"Using {args.threads} threads")
    logger.info(f"Using device: {args.device}")
    logger.info(f"OS: {platform.system()} {platform.machine()}")

    proteins = parse_fasta(args.input)

    if args.command == "nlr":
        proteins = hmmsearch(proteins, NLR_HMM_DB, threads=args.threads)

        if not args.retain:
            proteins = {k: v for k, v in proteins.items() if v.is_nlr()}

        nlrexpress(proteins, threads=args.threads)

        if args.coconat:
            predict_coils(proteins, args.device, args.batch_size, args.threads)

        for protein in proteins.values():
            protein.annotate_lrr(lrr_gap=args.lrr_gap, lrr_length=args.lrr_length)
            protein.merge_domains()
            protein.classify_nlr()

    elif args.command == "prr":
        proteins = hmmsearch(proteins, RLP_HMM_DB, threads=args.threads)
        tmbed(proteins, args.device, args.batch_size, args.threads)
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
