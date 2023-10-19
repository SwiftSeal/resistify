import subprocess
import sys
import logging
import os
from resistify.annotations import Annotation


def hmmsearch(input_file, sequences, temp_dir, data_dir):
    hmmsearch_db = os.path.join(data_dir, "nlrdb.hmm")
    output_file = os.path.join(temp_dir.name, "hmmsearch.out")

    cmd = [
        "hmmsearch",
        "--noali",
        "-E",
        "0.00001",
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
        logging.error(f"😞 Error running hmmsearch. Stderr of hmmsearch: {e.stderr}")
        sys.exit(1)
    except FileNotFoundError:
        logging.error(f"😞 hmmsearch not found. Have you installed it?")
        sys.exit(1)

    with open(output_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            line = line.split()
            sequence = line[0]
            domain = line[3]
            evalue = float(line[11])
            start = int(line[17])
            end = int(line[18])

            if start > end:
                continue

            if evalue > 0.00001:
                continue

            # MADA motif must be at the start of the sequence to be valid
            if domain == "MADA" and start != 1:
                continue

            if domain == "TIR_2":
                domain = "TIR"

            if domain == "Rx_N":
                domain = "CC"

            sequences[sequence].add_annotation(Annotation(domain, start, end))

    return sequences
