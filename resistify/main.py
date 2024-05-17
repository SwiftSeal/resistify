#!/usr/bin/env python3

import argparse
from rich_argparse import RichHelpFormatter
import tempfile
import sys
from .logging_setup import log, console
import os
from .annotations import *
from .hmmsearch import *
from .nlrexpress import *
from .utility import *

def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for identifying and classifying NLR resistance genes in plant genomes.
        """,
        formatter_class=RichHelpFormatter
    )
    parser.add_argument("-t", "--threads", help="Threads available to jackhmmer", default=2, type=int)
    parser.add_argument("--ultra", help="Run in ultra mode, non-NLRs will be retained", action="store_true")
    parser.add_argument("--chunksize", help="Number of sequences per split for jackhmmer", default=5, type=int)
    parser.add_argument("--evalue", help="E-value threshold for hmmsearch", default="0.00001")
    parser.add_argument("--lrr_gap", help="Minimum gap between LRR motifs", default=75, type=int)
    parser.add_argument("--lrr_length", help="Minimum number of LRR motifs to be considered an LRR domain", default=4, type=int)
    parser.add_argument("--duplicate_gap", help="Gap size (aa) to consider merging duplicate annotations", default=100, type=int)
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("outdir", help="Output directory")

    return parser.parse_args()

def main():
    log.info("Welcome to Resistify version 0.2.0!")

    args = parse_args()

    data_dir = os.path.join(os.path.dirname(__file__), "data")

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    temp_dir = tempfile.TemporaryDirectory()
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

        log.info(
            f"{len(classified_sequences)} sequences classified as potential NLRs!"
        )

    jackhmmer_input = save_fasta(
        classified_sequences, os.path.join(temp_dir.name, "jackhmmer_input.fa")
    )

    classified_sequences = jackhmmer(
        jackhmmer_input, classified_sequences, temp_dir, data_dir, args.chunksize, args.threads
    )

    # close the temporary directory
    temp_dir.cleanup()

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

    log.info("Thank you for using Resistify!")


if __name__ == "__main__":
    main()
