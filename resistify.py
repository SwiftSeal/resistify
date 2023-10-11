#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import logging
import os
from resistify.annotations import *
from resistify.hmmsearch import *
from resistify.nlrexpress import *

database_files = ["pfam.hmm", "superfamily.hmm", "smart.hmm", "gene3d.hmm", "cjid.hmm"]


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
            logging.error(f"😞 Database file {file}} does not exist")
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
    sequence_data = {}
    with open(path) as file:
        for record in SeqIO.parse(file, "fasta"):
            sequence_data[record.id] = record.seq
    return sequence_data


def save_fasta(sequence_data, path):
    with open(path, "w") as file:
        for sequence_id in sequence_data:
            sequence = sequence_data[sequence_id]
            file.write(f">{sequence_id}\n")
            file.write(f"{sequence.sequence}\n")
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
                    database_paths["superfamily"],
                    args.evalue,
                ),
                (input_fasta, temp_dir, database_path, "pfam", args.evalue),
                (input_fasta, temp_dir, database_path, "smart", args.evalue),
                (input_fasta, temp_dir, database_path, "cjid", args.evalue),
            ],
        )

    # fix accession names, merge and sort results into a single table
    results_file = save_fixed_accession(results, temp_dir, results_dir)

    sequence_annotations = parse_hmmer_table(results_file)

    for sequence_id in sequence_annotations:
        sequence = sequence_annotations[sequence_id]
        annotations = merge_and_sort(sequence.annotations)
        print(f"{sequence_id}\t{annotation_string(annotations)}")

    # TODO subset sequences based on annotations prior to motif prediction

    jackhmmer(input_fasta, temp_dir, "targetDB.fasta")

    os.path.join(temp_dir, "jackhmmer-1.hmm")

    jackhmmer_iteration_1 = parse_jackhmmer(
        os.path.join(temp_dir, "jackhmmer-1.hmm"), iteration=False
    )
    jackhmmer_iteration_2 = parse_jackhmmer(
        os.path.join(temp_dir, "jackhmmer-2.hmm"), iteration=True
    )

    input_data = generateInputFile(
        sequences, jackhmmer_iteration_1, jackhmmer_iteration_2
    )

    for predictor in motif_models.keys():
        result = predict_motif(sequences, input_data, predictor)

        result_index = 0
        for sequence in sequences:
            sequence_length = len(sequences[sequence])
            for i in range(len(sequences[sequence])):
                # make sure we are within the sequence bounds
                if i >= 5 and i < sequence_length - (motif_span_lengths[predictor] + 5):
                    value = round(result[result_index][1], 4)
                    if value > 0.8:
                        print(sequence, i, predictor, value)
                    result_index += 1


if __name__ == "__main__":
    main()
