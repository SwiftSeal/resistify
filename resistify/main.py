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
)
from resistify.hmmsearch import hmmsearch
from resistify.nlrexpress import nlrexpress
from resistify.coconat import coconat
from resistify.tmbed import tmbed

__version__ = "0.6.0"


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for identifying and classifying resistance genes in plant genomes.
        """,
        formatter_class=RichHelpFormatter,
    )

    # Add global arguments
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

    # Add subparsers
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
        "--batch",
        help="Number of sequences to process in parallel. This can help reduce memory usage.",
        default=None,
        type=int,
    )
    nlr_parser.add_argument(
        "--coconat",
        help="!EXPERIMENTAL! Path to the Coconat database. If provided, Coconat will be used to improve coiled-coil (CC) annotations.",
        action="store_true",
    )
    nlr_parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer. Default is 5.",
        default=5,
        type=int,
    )
    nlr_parser.add_argument(
        "--evalue",
        help="E-value threshold for hmmsearch. Default is 0.00001.",
        default="0.00001",
    )
    nlr_parser.add_argument(
        "--lrr_gap",
        help="Minimum gap (in amino acids) between LRR motifs. Default is 75.",
        default=75,
        type=int,
    )
    nlr_parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs required to be considered an LRR domain. Default is 4.",
        default=4,
        type=int,
    )
    nlr_parser.add_argument(
        "--duplicate_gap",
        help="Gap size (in amino acids) to consider merging duplicate annotations. Default is 100.",
        default=100,
        type=int,
    )
    nlr_parser.add_argument(
        "-o",
        "--outdir",
        help="Path to the output directory where results will be saved.",
        default=os.getcwd(),
    )
    nlr_parser.add_argument(
        "input",
        help="Path to the input FASTA file containing sequences to be analyzed.",
    )

    # PRR subparser
    prr_parser = subparsers.add_parser(
        "prr",
        help="Identify and classify PRR resistance genes.",
        formatter_class=RichHelpFormatter,
    )
    prr_parser.add_argument(
        "input",
        help="Path to the input FASTA file containing sequences to be analyzed.",
    )
    prr_parser.add_argument(
        "-o",
        "--outdir",
        help="Path to the output directory where results will be saved.",
        default=os.getcwd(),
    )
    prr_parser.add_argument(
        "--lrr_gap",
        help="Minimum gap (in amino acids) between LRR motifs. Default is 75.",
        default=75,
        type=int,
    )
    prr_parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs required to be considered an LRR domain. Default is 4.",
        default=4,
        type=int,
    )
    prr_parser.add_argument(
        "--duplicate_gap",
        help="Gap size (in amino acids) to consider merging duplicate annotations. Default is 100.",
        default=100,
        type=int,
    )
    prr_parser.add_argument(
        "--evalue",
        help="E-value threshold for hmmsearch. Default is 0.00001.",
        default="0.00001",
    )
    prr_parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer. Default is 5.",
        default=5,
        type=int,
    )
    # PRR-specific arguments can be added later

    return parser.parse_args()


def nlr(args, log):
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
        sequences = coconat(sequences)

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
    log.info("Searching for PRRs...")
    sequences = parse_fasta(args.input)
    sequences = hmmsearch(sequences, "prr", args.evalue)
    sequences = nlrexpress(
        sequences,
        "lrr",
        args.chunksize,
    )

    sequences = tmbed(
        sequences
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
