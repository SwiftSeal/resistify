#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import logging
import os
from resistify.annotations import *
from resistify.hmmsearch import *

database_paths = {
    "pfam": "resistify/data/pfam.hmm",
    "superfamily": "resistify/data/superfamily.hmm",
    "smart": "resistify/data/smart.hmm",
    "gene3d": "resistify/data/gene3d.hmm",
    "cjid": "resistify/data/abe3069_Data_S1.hmm",
}


def check_database_paths(database_paths):
    for source in database_paths:
        if os.path.exists(database_paths[source]):
            logging.info(
                f"😊 Database file for {source} found at {database_paths[source]}"
            )
        else:
            logging.error(
                f"😞 Database file for {source} does not exist at {database_paths[source]}"
            )
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

    check_database_paths(database_paths)

    results_dir = create_output_directory(args.outdir)

    with multiprocessing.Pool(5) as pool:
        results = pool.starmap(
            hmmsearch,
            [
                (args.fasta, "gene3d", database_paths["gene3d"], args.evalue),
                (
                    args.fasta,
                    "superfamily",
                    database_paths["superfamily"],
                    args.evalue,
                ),
                (args.fasta, "pfam", database_paths["pfam"], args.evalue),
                (args.fasta, "smart", database_paths["smart"], args.evalue),
                (args.fasta, "cjid", database_paths["cjid"], args.evalue),
            ],
        )

    results_file = save_fixed_accession(results, results_dir)

    sequences = parse_hmmer_table(results_file)

    for sequence_id in sequences:
        sequence = sequences[sequence_id]
        annotations = merge_and_sort(sequence.annotations)
        print(f"{sequence_id}\t{annotation_string(annotations)}")

    else:
        logging.error("😞 No command specified! Try resistify.py -h for help.")


if __name__ == "__main__":
    main()
