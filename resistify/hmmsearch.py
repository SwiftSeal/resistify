import subprocess
import sys
import logging
import os
from resistify.annotations import Annotation

def hmmsearch(input_file, sequences, temp_dir, data_dir, evalue):
    hmmsearch_db = os.path.join(data_dir, "nlrdb.hmm")
    output_file = os.path.join(temp_dir.name, "hmmsearch.out")

    cmd = [
        "hmmsearch",
        "--noali",
        "--domE",
        evalue,
        "--domtblout",
        output_file,
        hmmsearch_db,
        input_file,
    ]

    logging.info(f"😊 Running hmmsearch...")
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        logging.info(f"😊 hmmsearch completed successfully...")
    except subprocess.CalledProcessError as e:
        logging.error(f"😞 Error running hmmsearch:\nStderr: {e.stderr}\nStdout:{e.stdout}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error(f"😞 hmmsearch not found. Have you installed it?")
        sys.exit(1)

    sequences = parse_hmmsearch(output_file, sequences)

def parse_hmmsearch(output_file, sequences):

    with open(output_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.split()
            sequence = line[0]
            domain = line[3]
            evalue = float(line[11])
            score = float(line[13])
            start = int(line[17])
            end = int(line[18])

            if start > end:
                continue

            if domain == "TIR_2":
                domain = "TIR"

            if domain == "Rx_N" or domain == "cd14798":
                domain = "CC"

            # RPW8 being problematic - increase score threshold
            # what could possibly go wrong?
            if domain == "RPW8" and score < 20:
                continue

            sequences[sequence].add_annotation(Annotation(domain, start, end, evalue, score))

    return sequences
