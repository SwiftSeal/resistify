import subprocess
import sys
import os
import logging
from resistify.annotations import Annotation
from Bio import SearchIO
import tempfile

log = logging.getLogger(__name__)

def hmmsearch(sequences, evalue):
    hmmsearch_db = os.path.join(os.path.dirname(__file__), "data", "nlrdb.hmm")
    output_file = tempfile.NamedTemporaryFile()

    input_fasta = tempfile.NamedTemporaryFile()

    log.debug(f"Writing sequences to {input_fasta.name}")
    with open(input_fasta.name, "w") as f:
        for sequence in sequences:
            f.write(f">{sequence.id}\n{sequence.sequence}\n")
        
    cmd = [
        "hmmsearch",
        "--noali",
        "--domE",
        evalue,
        "--domtblout",
        output_file.name,
        hmmsearch_db,
        input_fasta.name,
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

    hmmsearch_results = {}
    for record in SearchIO.parse(output_file.name, "hmmsearch3-domtab"):
        for hit in record:
            for hsp in hit:
                # Necessary to avoid parsing errors
                if record.id == "TIR2":
                    record.id = "TIR"
                # Set a high threshold for RPW8 to avoid false positives
                if record.id == "RPW8" and hsp.bitscore < 20:
                    continue

                hmmsearch_results.setdefault(hit.id, []).append(
                    Annotation(
                        record.id, hsp.env_start, hsp.env_end, hsp.evalue, hsp.bitscore
                    )
                )
    
    for sequence in sequences:
        annotations = hmmsearch_results.get(sequence.id, [])
        for annotation in annotations:
            sequence.add_annotation(annotation)
    
    return sequences
