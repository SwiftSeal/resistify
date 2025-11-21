import os
import platform
import torch
import argparse
import logging
from pathlib import Path
from resistify.annotation import save_results
from resistify.hmmer import hmmsearch
from resistify.nlrexpress import nlrexpress
from resistify.coconat import predict_coils

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


def parse_args():
    parser = argparse.ArgumentParser(
        description="Resistify: A tool for NLR annotation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", type=Path, help="Path to the input FASTA file")
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="Output directory for results.",
    )
    parser.add_argument(
        "--coconat",
        action="store_true",
        help="Enable CoCoNat coiled-coil (CC) domain predictions. CCs will be predicted in the N-terminal region only.",
    )
    parser.add_argument(
        "--device",
        type=str,
        default=DEFAULT_DEVICE,
        help="Torch device to be used - can be 'cpu', 'cuda', or 'mps'.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=8,
        help="Batch size for ESM predictions - reduce this value for lower memory usage.",
    )
    parser.add_argument(
        "--pfam",
        type=Path,
        default=Path(os.path.dirname(__file__)) / "data" / "db.hmm",
        help="Path to the Pfam database file - by default a smaller NLR-specific database is used.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")

    logger.info("Resistify 2.0.0")
    logger.info(f"Python version: {platform.python_version()}")
    logger.info(f"PyTorch version: {torch.__version__}")
    logger.info(f"Using device: {args.device}")
    logger.info(f"OS: {platform.system()} {platform.machine()}")

    proteins = hmmsearch(args.fasta, args.pfam)

    nlrexpress(proteins)

    if args.coconat:
        predict_coils(proteins, args.device)

    for protein in proteins.values():
        protein.annotate_lrr()
        protein.merge_domains()
        protein.classify_nlr()

    save_results(proteins, args.output_dir)

    logger.info("Thank you for using Resistify!")
    logger.info("https://doi.org/10.1177/11779322241308944")
