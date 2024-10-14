import subprocess
import sys
import os
import logging
from resistify.annotations import Annotation
from resistify.utility import split_fasta
from Bio import SearchIO

log = logging.getLogger(__name__)

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

    sequences = parse_hmmsearch(output_file, sequences)
    
    return sequences


def parse_hmmsearch(output_file, sequences):
    for record in SearchIO.parse(output_file, "hmmsearch3-domtab"):
        for hit in record:
            for hsp in hit:
                # Necessary to avoid parsing errors
                if record.id == "TIR2":
                    record.id = "TIR"
                # Set a high threshold for RPW8 to avoid false positives
                if record.id == "RPW8" and hsp.bitscore < 20:
                    continue
                sequences[hit.id].add_annotation(
                    Annotation(
                        record.id, hsp.env_start, hsp.env_end, hsp.evalue, hsp.bitscore
                    )
                )
    return sequences
