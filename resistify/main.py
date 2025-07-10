import argparse
from pathlib import Path
from resistify.__version__ import __version__
from resistify._loguru import logger
from resistify.utility import get_threads
from resistify.nlr import nlr
from resistify.prr import prr
from resistify.draw import draw
from resistify.download_models import download_models
from resistify.check_dependencies import check_dependencies

DIAGRAM = """
e.g., resistify nlr proteins.fa -o results/───────────────────┐
                 │       │                    ┌───────────────▼────────────────┐
    ┌────────────┼───────┘                    │* NLR identification            │
    │            │    ┌────────────────────┐  │* NLR classification            │
    │          ┌─▼─┐  │* hmmsearch         │  │* CC, TIR, RPW8                 │
    │       ┌─►│NLR├─►│* NLRexpress        ├─►│* NB-ARC domains                │
    │       │  └───┘  │* CoCoNat (optional)│  │* Motif and domain annotations  │
┌───▼─────┐ │         └────────────────────┘  └────────────────────────────────┘
│Resistify├─┤
└─────────┘ │         ┌────────────────────┐  ┌────────────────────────────────┐
            │  ┌───┐  │* hmmsearch         │  │* RLK/RLP identification        │
            └─►│PRR├─►│* NLRexpress (LRR)  ├─►│* Extracellular classification  │
               └───┘  │* TMbed             │  │* Signal peptide identification │
                      └────────────────────┘  │* Motif and domain annotations  │
                                              └────────────────────────────────┘
"""


def main():
    parser = argparse.ArgumentParser(
        description="Resistify: A tool for identifying and classifying plant resistance genes.",
        epilog=DIAGRAM,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Resistify v{__version__}",
        help="Show the version of Resistify.",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Enable verbose logging."
    )

    subparsers = parser.add_subparsers(
        dest="subcommand", title="subcommands", description="Available subcommands"
    )

    nlr_parser = subparsers.add_parser(
        "nlr",
        help="Identify and classify NLR sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    nlr_parser.add_argument(
        "infile",
        type=Path,
        default=Path("."),
        help="Input FASTA file containing sequences to analyse.",
    )
    nlr_parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("."),
        help="Directory to save results.",
    )
    nlr_parser.add_argument(
        "--coconat",
        action="store_true",
        help="Run CoCoNat to identify additional CC domains.",
    )
    nlr_parser.add_argument(
        "--retain",
        action="store_true",
        help="Retain all sequences, not just those with NB-ARC domains.",
    )
    nlr_parser.add_argument(
        "--lrr-gap",
        type=int,
        default=75,
        help="Minimum gap (in amino acids) between LRR motifs.",
    )
    nlr_parser.add_argument(
        "--lrr-length",
        type=int,
        default=4,
        help="Minimum number of LRR motifs required to be considered an LRR domain.",
    )
    nlr_parser.add_argument(
        "--duplicate-gap",
        type=int,
        default=100,
        help="Gap size (in amino acids) to consider merging duplicate annotations.",
    )
    nlr_parser.add_argument(
        "--chunksize",
        type=int,
        default=None,
        help="Number of sequences per split for NLRexpress.",
    )
    nlr_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=get_threads(),
        help="Number of threads available for NLRexpress.",
    )

    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR sequences.",
        description="This command will identify and classify pattern recognition receptor (PRR) sequences from a given FASTA file. First, TMbed is used to predict transmembrane domains which is used as evidence for a receptor-like architecture. Then, HMM profiles of common extracellular domains and NLRexpress are used to annotate the extracellular domains. Finally, the sequences are classified into RLPs and RLKs based on the presence of an internal kinase domain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    prr_parser.add_argument(
        "infile", type=Path, help="Input FASTA file containing sequences to analyze."
    )
    prr_parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("."),
        help="Directory to save results.",
    )
    prr_parser.add_argument(
        "--lrr-gap",
        type=int,
        default=75,
        help="Minimum gap (in amino acids) between LRR motifs.",
    )
    prr_parser.add_argument(
        "--lrr-length",
        type=int,
        default=4,
        help="Minimum number of LRR motifs required to be considered an LRR domain.",
    )
    prr_parser.add_argument(
        "--duplicate-gap",
        type=int,
        default=100,
        help="Gap size (in amino acids) to consider merging duplicate annotations.",
    )
    prr_parser.add_argument(
        "--chunksize",
        type=int,
        default=5,
        help="Number of sequences per split for NLRexpress.",
    )
    prr_parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=get_threads(),
        help="Number of threads available for NLRexpress.",
    )

    draw_parser = subparsers.add_parser(
        "draw",
        help="Draw NLRs or PRRs from results directory.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    draw_parser.add_argument(
        "results_dir", type=Path, help="Directory containing results files."
    )
    draw_parser.add_argument(
        "-o",
        "--outfile",
        type=Path,
        default=Path("resistify_plot.png"),
        help="Output file for the plot.",
    )
    draw_parser.add_argument(
        "--query",
        type=str,
        default=None,
        help="Specific sequence to highlight in the plot.",
    )
    draw_parser.add_argument(
        "--width", type=int, default=10, help="Width of the plot in inches."
    )
    draw_parser.add_argument(
        "--height", type=int, default=5, help="Height of the plot in inches."
    )
    draw_parser.add_argument(
        "--hide-motifs", action="store_true", help="Hide motifs in the plot."
    )
    draw_parser.add_argument(
        "--dpi", type=int, default=300, help="DPI for the output image."
    )

    subparsers.add_parser("download_models", help="Download models manually.")

    subparsers.add_parser(
        "check_dependencies",
        help="Check if all dependencies are installed.",
    )

    args = parser.parse_args()

    if args.verbose:
        logger.update_level("DEBUG")

    logger.info(f"Welcome to Resistify v{__version__}!")
    logger.info("Need help? Visit https://github.com/SwiftSeal/Resistify")

    if args.subcommand == "nlr":
        nlr(
            infile=args.infile,
            outdir=args.outdir,
            run_coconat=args.coconat,
            retain=args.retain,
            lrr_gap=args.lrr_gap,
            lrr_length=args.lrr_length,
            duplicate_gap=args.duplicate_gap,
            chunksize=args.chunksize,
            threads=args.threads,
        )
    elif args.subcommand == "prr":
        prr(
            infile=args.infile,
            outdir=args.outdir,
            lrr_gap=args.lrr_gap,
            lrr_length=args.lrr_length,
            duplicate_gap=args.duplicate_gap,
            chunksize=args.chunksize,
            threads=args.threads,
        )
    elif args.subcommand == "draw":
        draw(
            results_dir=args.results_dir,
            outfile=args.outfile,
            query=args.query,
            width=args.width,
            height=args.height,
            hide_motifs=args.hide_motifs,
            dpi=args.dpi,
        )
        # Exit here, no need for printing help
        return
    elif args.subcommand == "download_models":
        download_models()
    elif args.subcommand == "check_dependencies":
        check_dependencies()
        return

    logger.info("Thank you for using Resistify!")
    logger.info("If you used Resistify in your research, please cite the following:")
    logger.info(" - Resistify: https://doi.org/10.1177/11779322241308944")
    logger.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    if args.subcommand == "prr":
        logger.info(" - TMbed: https://doi.org/10.1186/s12859-022-04873-x")
    elif args.subcommand == "nlr" and args.coconat:
        logger.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495")


if __name__ == "__main__":
    main()
