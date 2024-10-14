import argparse
from rich_argparse import RichHelpFormatter
import logging
from rich.logging import RichHandler
import sys
import os
from resistify.utility import (
    prepare_temp_directory,
    create_output_directory,
    parse_fasta,
    save_fasta,
    result_table,
    domain_table,
    motif_table,
    extract_nbarc,
)
from resistify.hmmsearch import hmmsearch
from resistify.nlrexpress import jackhmmer, motif_models, predict_motif
from resistify.annotations import Sequence

__version__ = "0.4.0"

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for identifying and classifying NLR resistance genes in plant genomes.
        """,
        formatter_class=RichHelpFormatter,
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    parser.add_argument(
        "-t", "--threads", help="Threads available to jackhmmer", default=2, type=int
    )
    parser.add_argument("--debug", help="Enable debug logging", action="store_true")
    parser.add_argument(
        "--ultra",
        help="Run in ultra mode, non-NLRs will be retained",
        action="store_true",
    )
    parser.add_argument(
        "--chunksize",
        help="Number of sequences per split for jackhmmer",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--evalue", help="E-value threshold for hmmsearch", default="0.00001"
    )
    parser.add_argument(
        "--lrr_gap", help="Minimum gap between LRR motifs", default=75, type=int
    )
    parser.add_argument(
        "--lrr_length",
        help="Minimum number of LRR motifs to be considered an LRR domain",
        default=4,
        type=int,
    )
    parser.add_argument(
        "--duplicate_gap",
        help="Gap size (aa) to consider merging duplicate annotations",
        default=100,
        type=int,
    )
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("outdir", help="Output directory")

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

    data_dir = os.path.join(os.path.dirname(__file__), "data")
    temp_dir = prepare_temp_directory(data_dir)
    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    log.debug(f"Temporary directory created at {temp_dir.name}")

    # save the fasta with stripped headers
    hmmer_input = save_fasta(sequences, os.path.join(temp_dir.name, "hmmer_input.fa"))

    sequences = hmmsearch(hmmer_input, sequences, temp_dir, data_dir, args.evalue)

    for sequence in sequences:
        sequences[sequence].merge_annotations()
        sequences[sequence].classify()

    # subset sequences based on classification
    if args.ultra:
        # do not subset sequences based on classification
        log.info(f"Running in ultra mode!")
        classified_sequences = sequences
    else:
        # subset sequences based on classification
        classified_sequences = {
            sequence: sequences[sequence]
            for sequence in sequences
            if sequences[sequence].classification is not None
        }

        if not classified_sequences:
            log.info(f"No sequences classified as potential NLRs!")
            sys.exit(0)

        log.info(f"{len(classified_sequences)} sequences classified as potential NLRs!")

    jackhmmer_input = save_fasta(
        classified_sequences, os.path.join(temp_dir.name, "jackhmmer_input.fa")
    )

    classified_sequences = jackhmmer(
        jackhmmer_input,
        classified_sequences,
        temp_dir,
        data_dir,
        args.chunksize,
        args.threads,
    )

    # close the temporary directory
    # temp_dir.cleanup()

    # predict and add motifs to sequences
    # perhaps move all of this into a function rather than iterating
    for predictor in motif_models.keys():
        predict_motif(classified_sequences, predictor, data_dir)

    for sequence in classified_sequences:
        classified_sequences[sequence].reclassify(args.lrr_gap, args.lrr_length)

    log.info(f"Saving results to {results_dir}...")
    result_table(classified_sequences, results_dir)
    domain_table(classified_sequences, results_dir)
    motif_table(classified_sequences, results_dir)
    extract_nbarc(classified_sequences, results_dir)
    save_fasta(
        classified_sequences, os.path.join(results_dir, "nlr.fasta"), nlr_only=True
    )

    log.info("Thank you for using Resistify!")


if __name__ == "__main__":
    main()
