#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import logging
import os
from resistify.annotations import *
from resistify.hmmsearch import *
from resistify.nlrexpress import *
from resistify.utility import *


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

    hmmer_db = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data/nlrdb.hmm")

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    temp_dir = tempfile.TemporaryDirectory()

    # save the fasta with stripped headers
    input_fasta = save_fasta(sequences, os.path.join(temp_dir.name, "input.fasta"))

    sequences = hmmsearch(input_fasta, sequences, temp_dir, hmmer_db)

    for sequence in sequences:
        sequences[sequence].merge_annotations()
        sequences[sequence].classify()


    """
    jackhmmer(input_fasta, temp_dir, database_path)

    jackhmmer_iteration_1 = parse_jackhmmer(
        os.path.join(temp_dir.name, "jackhmmer-1.hmm"), iteration=False
    )
    jackhmmer_iteration_2 = parse_jackhmmer(
        os.path.join(temp_dir.name, "jackhmmer-2.hmm"), iteration=True
    )

    sequences = prepare_jackhmmer_data(
        sequences, jackhmmer_iteration_1, jackhmmer_iteration_2
    )

    # close the temporary directory
    temp_dir.cleanup()

    # predict and add motifs to sequences
    # perhaps move all of this into a function rather than iterating
    for predictor in motif_models.keys():
        predict_motif(sequences, predictor)

    # print a table of the result
    for sequence in sequences:
        sequence_string = sequences[sequence].annotation_string()
        downstream_lrrs = len(sequences[sequence].downstream_lrr())
        print(sequence, sequence_string, downstream_lrrs)
    """


if __name__ == "__main__":
    main()
