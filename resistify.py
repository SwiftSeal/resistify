import argparse
import sys
import logging
from resistify.annotations import Sequence, Annotation, parse_hmmer_table, classifications
from resistify.hmmsearch import hmmsearch

database_paths = {
    "pfam": "data/pfam.hmm",
    "superfamily": "data/superfamily.hmm",
    "smart": "data/smartsmart.hmm",
    "gene3d": "data/gene3d.hmm"
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse HMMER output to generate a table of annotations."
    )
    parser.add_argument("--verbose", help="Verbose output", action="store_true")
    subparsers = parser.add_subparsers(help = "sub-command help", dest = "command")
    search = subparsers.add_parser("search", help = "Search a protein sequence against a database of HMMs")
    search.add_argument("fasta", help = "Protein sequences to search")
    #search.add_argument("--db", help = " Path to HMM databases", type = str, default = "./resistify/data")
    search.add_argument("-e", "--evalue", help = "E-value threshold", type = float, default = 1e-5)

    annotate = subparsers.add_parser("annotate", help = "Generate a table of annotations from HMMER output")
    annotate.add_argument("input", nargs = "+", type = str, help = "HMMER output files")
    annotate.add_argument("-e", "--evalue", help = "E-value threshold", type = float, default = 1e-5)

    return parser.parse_args()

def main():
    args = parse_args()

    if args.verbose:
        logging.basicConfig(level = logging.INFO, stream=sys.stderr)
    else:
        logging.basicConfig(level = logging.WARNING, stream=sys.stderr)

    if args.command == "search":
        for file in args.fasta:
            hmmsearch(file, "pfam", database_paths["pfam"], args.evalue)
            hmmsearch(file, "superfamily", database_paths["superfamily"], args.evalue)
            hmmsearch(file, "smart", database_paths["smart"], args.evalue)
            hmmsearch(file, "gene3d", database_paths["gene3d"], args.evalue)
    elif args.command == "annotate":
        sequences = {}
        for file in args.input:
            sequences = parse_hmmer_table(file, sequences, args.evalue)
        for sequence_name in sequences:
            sequence = sequences[sequence_name]
            print(f"{sequence.name}\t{sequence.length}\t{sequence.annotation_string()}")
    else:
        logging.error("No command specified! Try resistify.py -h for help.")
    """
    sequences = parse_hmmer_table(args.input, args.evalue)

    for sequence_name in sequences:
        sequence = sequences[sequence_name]
        print(f"{sequence.name}\t{sequence.length}\t{sequence.annotation_string()}")
    """

if __name__ == "__main__":
    main()