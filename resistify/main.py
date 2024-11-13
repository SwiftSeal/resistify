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
from resistify.annotations import Sequence
from resistify.coconat import coconat
from resistify.tmbed import tmbed

__version__ = "0.5.2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for identifying and classifying NLR resistance genes in plant genomes.
        """,
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version number and exit.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        help="Total number of threads available for jackhmmer multiprocessing. By default, all available threads will be used.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--debug", help="Enable debug logging for detailed output.", action="store_true"
    )
    parser.add_argument(
        "--ultra",
        help="Run in ultra mode to retain non-NLR sequences.",
        action="store_true",
    )
    parser.add_argument(
        "--batch",
        help="Number of sequences to process in parallel. This can help reduce memory usage.",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--coconat",
        help="!EXPERIMENTAL! Path to the Coconat database. If provided, Coconat will be used to improve coiled-coil (CC) annotations.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer. Default is 5.",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--evalue",
        help="E-value threshold for hmmsearch. Default is 0.00001.",
        default="0.00001",
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
        "input",
        help="Path to the input FASTA file containing sequences to be analyzed.",
    )
    parser.add_argument(
        "outdir", help="Path to the output directory where results will be saved."
    )

    return parser.parse_args()


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

    if args.coconat:
        log.info(
            f"CoCoNat database provided - this will be used to improve CC annotations."
        )

    # Calculate threads to use
    if args.threads is None:
        # Use all available threads by default
        thread_count = len(os.sched_getaffinity(0))
        log.debug(f"Using {thread_count} threads by default.")
    else:
        thread_count = args.threads

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    sequences = hmmsearch(sequences, args.evalue)

    for sequence in sequences:
        sequence.classify()

    # subset sequences based on classification
    if args.ultra:
        # do not subset sequences based on classification
        log.info(f"Running in ultra mode!")
    else:
        # subset sequences based on classification
        sequences = [
            sequence for sequence in sequences if sequence.classification is not None
        ]

        if not sequences:
            log.info(f"No sequences classified as potential NLRs!")
            sys.exit(0)

        log.info(f"{len(sequences)} sequences classified as potential NLRs!")

    if args.batch is None:
        batch_size = len(sequences)
    else:
        batch_size = args.batch

    batches = [
        sequences[i : i + batch_size]
        for i in range(0, len(sequences), batch_size)
    ]

    for batch in batches:
        log.info(f"Processing batch of {len(batch)} sequences...")
        if args.coconat:
            log.info(f"Running CoCoNat...")
            batch = coconat(batch, args.coconat)
            for sequence in batch:
                sequence.identify_cc_domains()
                sequence.classify()
    
        batch = nlrexpress(
            batch,
            args.chunksize,
            thread_count,
        )
    
    sequences = [sequence for batch in batches for sequence in batch]

    for sequence in sequences:
        sequence.identify_lrr_domains(args.lrr_gap, args.lrr_length)
        sequence.classify()
        sequence.merge_annotations(args.duplicate_gap)
    
    sequences = tmbed(sequences)

    log.info(f"Saving results to {results_dir}...")
    result_table(sequences, results_dir)
    annotation_table(sequences, results_dir)
    domain_table(sequences, results_dir)
    motif_table(sequences, results_dir)
    extract_nbarc(sequences, results_dir)
    save_fasta(
        sequences, os.path.join(results_dir, "nlr.fasta"), nlr_only=True
    )
    if args.coconat:
        coconat_table(sequences, results_dir)

    log.info("Thank you for using Resistify!")
    log.info("If you used Resistify in your research, please cite the following:")
    log.info(" - Resistify: https://doi.org/10.1101/2024.02.14.580321")
    log.info(" - NLRexpress: https://doi.org/10.3389/fpls.2022.975888")
    log.info(" - CoCoNat: https://doi.org/10.1093/bioinformatics/btad495 (if used)")


if __name__ == "__main__":
    main()
