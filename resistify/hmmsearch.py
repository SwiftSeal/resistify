import subprocess
import sys
import os
import logging
from resistify.annotations import Annotation
from resistify.utility import split_fasta
from Bio import SearchIO

log = logging.getLogger(__name__)

def hmmsearch_subprocess(input, output, db, evalue):
    cmd = [
        "hmmsearch",
        "--noali",
        "--domE",
        evalue,
        "--domtblout",
        output,
        db,
        input,
    ]

    log.info(f"Running hmmsearch...")
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        log.info(f"hmmsearch completed successfully...")
    except subprocess.CalledProcessError as e:
        log.error(f"Error running hmmsearch:\nStderr: {e.stderr}\nStdout:{e.stdout}")
        sys.exit(1)
    except FileNotFoundError:
        log.error(f"hmmsearch not found. Have you installed it?")
        sys.exit(1)


def hmmsearch(input_file, sequences, temp_dir, data_dir, evalue):
    hmmsearch_db = os.path.join(data_dir, "nlrdb.hmm")
    output_file = os.path.join(temp_dir.name, "hmmsearch.out")

    if len(sequences) >= 100000:
        log.debug(f"Splitting input file into chunks")
        fastas = split_fasta(input_file, 50000, temp_dir)
        for fasta in fastas:
            hmmsearch_subprocess(fasta, f"{fasta}.out", hmmsearch_db, evalue)
            sequences = parse_hmmsearch(f"{fasta}.out", sequences)
    else:
        hmmsearch_subprocess(input_file, output_file, hmmsearch_db, evalue)
        sequences = parse_hmmsearch(output_file, sequences)
    return sequences


def parse_hmmsearch(output_file, sequences):
    for record in SearchIO.parse(output_file, "hmmsearch3-domtab"):
        for hit in record:
            for hsp in hit:
                # Set a high threshold for RPW8 to avoid false positives
                if record.id == "RPW8" and hsp.bitscore < 20:
                    continue
                sequences[hit.id].add_annotation(
                    Annotation(
                        record.id, hsp.env_start, hsp.env_end, hsp.evalue, hsp.bitscore
                    )
                )
    return sequences
