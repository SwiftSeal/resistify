#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import logging
import os
from resistify.annotations import *
from resistify.hmmsearch import *
from resistify.nlrexpress import *

database_files = ["pfam.hmm", "superfamily.hmm", "smart.hmm", "gene3d.hmm", "gene3d.tsv", "cjid.hmm"]


def check_database(database_path):
    # check if database_path exists
    if os.path.exists(database_path):
        logging.info(f"😊 Database directory exists at {database_path}")
    else:
        logging.error(f"😞 Database directory does not exist at {database_path}")
        sys.exit(1)

    for file in database_files:
        if os.path.exists(os.path.join(database_path, file)):
            logging.info(f"😊 Database file {file} exists")
        else:
            logging.error(f"😞 Database file {file} does not exist")
            sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Resistify is a tool for predicting resistance genes in plant genomes.
        Domains are predicted using hmmsearch against multiple databases which are classified as common resistance gene related domains.
        Sequences are then classified based on the presence and order of these domains.
        """
    )
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    parser.add_argument("--evalue", help="E-value threshold", default="0.00001")
    parser.add_argument("--database_path", help="Path to the database directory")
    parser.add_argument("input", help="Input FASTA file")
    parser.add_argument("outdir", help="Output directory")

    return parser.parse_args()


def create_output_directory(outdir):
    try:
        expanded_outdir = os.path.expanduser(os.path.expandvars(outdir))
        os.makedirs(expanded_outdir, exist_ok=True)
        logging.info(f"😊 Output directory created at {expanded_outdir}")
        return expanded_outdir
    except OSError as e:
        logging.error(f"😞 Error creating output directory: {e}")
        sys.exit(1)


def parse_fasta(path):
    sequences = {}
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences[record.id] = Sequence(record.id, record.seq)
    return sequences


def save_fasta(sequences, path):
    with open(path, "w") as file:
        for sequence in sequences:
            file.write(f">{sequence}\n")
            file.write(f"{sequences[sequence].sequence}\n")
    return path


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

    if args.database_path:
        database_path = args.database_path
    else:
        # default to the database directory in the same directory as resistify.py
        database_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "data"
        )

    check_database(database_path)

    results_dir = create_output_directory(args.outdir)

    sequences = parse_fasta(args.input)

    # create a temporary directory to store tempfiles
    temp_dir = tempfile.TemporaryDirectory()

    # save the fasta with stripped headers
    input_fasta = save_fasta(sequences, os.path.join(temp_dir.name, "input.fasta"))

    with multiprocessing.Pool(5) as pool:
        results = pool.starmap(
            hmmsearch,
            [
                (input_fasta, temp_dir, database_path, "gene3d", args.evalue),
                (
                    input_fasta,
                    temp_dir,
                    database_path,
                    "superfamily",
                    args.evalue,
                ),
                (input_fasta, temp_dir, database_path, "pfam", args.evalue),
                (input_fasta, temp_dir, database_path, "smart", args.evalue),
                (input_fasta, temp_dir, database_path, "cjid", args.evalue),
            ],
        )

    # fix accession names, merge and sort results into a single table
    results_file = save_fixed_accession(results, temp_dir, database_path, results_dir)

    sequences = parse_hmmer_table(sequences, results_file)

    for sequence in sequences:
        sequences[sequence].merge_annotations()
        print(sequence, sequences[sequence].annotation_string())

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
        lrrs = sequences[sequence].motifs["LxxLxL"]
        print(sequence, sequence_string, len(lrrs))

    


if __name__ == "__main__":
    main()
