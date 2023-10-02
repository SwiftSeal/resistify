import argparse
from resistify.annotations import Annotation, Sequence, classifications, parse_hmmer_table

def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse HMMER output to generate a table of annotations."
    )
    parser.add_argument("-i", "--input", help="HMMER output file")
    parser.add_argument(
        "-e", "--evalue", help="E-value threshold", type=float, default=1e-2
    )
    return parser.parse_args()


def main():
    args = parse_args()

    sequences = parse_hmmer_table(args.input, args.evalue)

    for sequence_name in sequences:
        sequence = sequences[sequence_name]
        print(f"{sequence.name}\t{sequence.length}\t{sequence.annotation_string()}")

if __name__ == "__main__":
    main()