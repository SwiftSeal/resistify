import os
import argparse
from argparse import RawDescriptionHelpFormatter
from resistify.__version__ import __version__

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


def validate_input_file(filepath: str) -> str:
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
        formatter_class=RawDescriptionHelpFormatter,
        description="A tool for identifying and classifying resistance genes in plant genomes.",
        epilog=DIAGRAM,
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
        help="If enabled, CoCoNat will be used to improve coiled-coil (CC) annotations.",
        action="store_true",
    )
    add_common_args(nlr_parser)

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes.",
    )
    add_common_args(prr_parser)

    # Download models subparser
    subparsers.add_parser(
        "download_models",
        help="Download the models required for CoCoNat and TMbed. These will be stored in the default $HF_HOME and $TORCH_HOME directories. Otherwise, these will be downloaded automatically when required.",
    )

    # Draw subparser
    draw_parser = subparsers.add_parser(
        "draw",
        help="Draw domain structure for target gene(s) from a given results directory.",
    )
    draw_parser.add_argument(
        "--query",
        help="Comma-separated list of sequence names to plot. If not provided, all sequences will be plotted.",
        default=None,
    )
    draw_parser.add_argument(
        "results_dir",
        help="Path to the results directory.",
    )
    draw_parser.add_argument(
        "-o",
        "--output",
        help="Path to the output plot. Extension should be .png, .pdf, or .svg. Default is domain_plot.png",
        default="domain_plot.png",
        type=str,
    )
    draw_parser.add_argument(
        "--width",
        help="Width of the output plot in inches.",
        type=float,
        default=12,
    )
    draw_parser.add_argument(
        "--height",
        help="Height of the output plot in inches.",
        type=float,
        default=6,
    )
    draw_parser.add_argument(
        "--hide-motifs",
        help="Hide motifs in the output plot.",
        action="store_true",
    )

    return parser.parse_args(args)
