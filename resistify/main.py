#!/usr/bin/env python3

import argparse
import tempfile
import sys
import logging
import os
from .annotations import *
from .hmmsearch import *
from .nlrexpress import *
from .utility import *


def main():
    args = parse_args()

    if args.verbose:
        logging.basicConfig(
            level=logging.DEBUG, stream=sys.stderr, format="%(asctime)s - %(message)s"
        )
    else:
        logging.basicConfig(
            level=logging.INFO, stream=sys.stderr, format="%(asctime)s - %(message)s"
        )

    data_dir = os.path.join(os.path.dirname(__file__), "data")

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    temp_dir = tempfile.TemporaryDirectory()

    # save the fasta with stripped headers
    hmmer_input = save_fasta(sequences, os.path.join(temp_dir.name, "hmmer_input.fa"))

    sequences = hmmsearch(hmmer_input, sequences, temp_dir, data_dir)

    for sequence in sequences:
        sequences[sequence].merge_annotations()
        sequences[sequence].classify()

    # subset sequences based on classification
    classified_sequences = {
        sequence: sequences[sequence]
        for sequence in sequences
        if sequences[sequence].classification is not None
    }

    if not classified_sequences:
        logging.info(f"😞 No sequences classified as potential NLRs!")
        sys.exit(0)

    logging.info(
        f"😊 {len(classified_sequences)} sequences classified as potential NLRs!"
    )

    jackhmmer_input = save_fasta(
        classified_sequences, os.path.join(temp_dir.name, "jackhmmer_input.fa")
    )

    classified_sequences = jackhmmer(
        jackhmmer_input, classified_sequences, temp_dir, data_dir
    )

    # close the temporary directory
    temp_dir.cleanup()

    # predict and add motifs to sequences
    # perhaps move all of this into a function rather than iterating
    for predictor in motif_models.keys():
        predict_motif(classified_sequences, predictor, data_dir)

    for sequence in classified_sequences:
        classified_sequences[sequence].reclassify()

    logging.info(f"😊 Saving results to {results_dir}...")
    result_table(classified_sequences, results_dir)
    domain_table(classified_sequences, results_dir)
    motif_table(classified_sequences, results_dir)
    extract_nbarc(classified_sequences, results_dir)
    logging.info(f"😊 Thank you for using Resistify!")


if __name__ == "__main__":
    main()
