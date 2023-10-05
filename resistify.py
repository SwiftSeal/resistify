#!/usr/bin/env python3

import argparse
import multiprocessing
import sys
import logging
import os
from resistify.annotations import (
    Sequence,
    Annotation,
    parse_hmmer_table,
    classifications,
)
from resistify.hmmsearch import hmmsearch, print_fixed_accession

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
        description="Parse HMMER output to generate a table of annotations."
    )
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    subparsers = parser.add_subparsers(help="sub-command help", dest="command")
    search = subparsers.add_parser(
        "search", help="Search a protein sequence against a database of HMMs"
    )
    search.add_argument("fasta", help="Protein sequences to search")
    search.add_argument(
        "-e", "--evalue", help="E-value threshold", type=str, default="0.00001"
    )

    annotate = subparsers.add_parser(
        "annotate", help="Generate a table of annotations from HMMER output"
    )
    annotate.add_argument("input", nargs="+", type=str, help="HMMER output files")
    annotate.add_argument(
        "-e", "--evalue", help="E-value threshold", type=str, default="0.00001"
    )

    return parser.parse_args()

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

    if args.command == "search":
        check_database_paths(database_paths)

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

        for result in results:
            print_fixed_accession(result)


    elif args.command == "annotate":
        sequences = {}
        for file in args.input:
            sequences = parse_hmmer_table(file, sequences, args.evalue)
        for sequence_name in sequences:
            sequence = sequences[sequence_name]
            print(f"{sequence.name}\t{sequence.length}\t{sequence.annotation_string()}")
    else:
        logging.error("😞 No command specified! Try resistify.py -h for help.")


if __name__ == "__main__":
    main()
