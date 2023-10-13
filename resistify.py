#!/usr/bin/env python3

import argparse
import tempfile
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
    jackhmmer_db = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data/nlrexpress.fasta")

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    temp_dir = tempfile.TemporaryDirectory()

    # save the fasta with stripped headers
    hmmer_input = save_fasta(sequences, os.path.join(temp_dir.name, "hmmer_input.fa"))

    sequences = hmmsearch(hmmer_input, sequences, temp_dir, hmmer_db)

    for sequence in sequences:
        sequences[sequence].merge_annotations()
        sequences[sequence].classify()
        print(sequence, sequences[sequence].classification)

    # save only sequences with a classification
    classified = [sequence for sequence in sequences if sequences[sequence].classification is not None]
    jackhmmer_input = save_fasta(sequences, os.path.join(temp_dir.name, "jackhmmer_input.fa"), classified)
    
    sequence = jackhmmer(jackhmmer_input, sequences, temp_dir, jackhmmer_db)

    # close the temporary directory
    temp_dir.cleanup()

    # predict and add motifs to sequences
    # perhaps move all of this into a function rather than iterating
    for predictor in motif_models.keys():
        predict_motif(sequences, predictor)


if __name__ == "__main__":
    main()
