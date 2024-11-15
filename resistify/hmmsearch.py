import subprocess
import sys
import os
import logging
from resistify.annotations import Annotation
from Bio import SearchIO
import tempfile

log = logging.getLogger(__name__)

accession_families = {
    "PF00931": "NB-ARC",
    "PF01582": "TIR",
    "PF05659": "RPW8",
    "PF13676": "TIR",
    "PF18052": "CC",
    "PF20160": "C-JID",
    "MADA": "MADA",
    "PF00069": "PKinase",
    "PF01453": "G-LecRLK",
    "PF00954": "G-LecRLK",
    "PF08276": "G-LecRLK",
    "PF00024": "G-LecRLK",
    "PF14295": "G-LecRLK",
    "PF00139": "L-LecRLK",
    "PF00059": "C-LecRLK",
    "PF13947": "WAK",
    "PF14380": "WAK",
    "PF00008": "WAK",
    "PF08488": "WAK",
    "PF07645": "WAK",
    "PF12662": "WAK",
    "PF12947": "WAK",
    "PF11721": "CrRLK1L",
    "PF12819": "CrRLK1L",
    "PF01476": "LysM",
    "PF01657": "CRK",
    "PF00314": "Thaumatin",
    "PF13540": "CR-like",
    "PF19160": "SPARK",
    "PF00704": "GH18",
    "PF00182": "GH19",
    "PF00188": "CAP",
    "PF16101": "PRIMA1"
}


def hmmsearch(sequences, evalue):
    hmmsearch_db = os.path.join(os.path.dirname(__file__), "data", "nlrdb.hmm")
    output_file = tempfile.NamedTemporaryFile()

    input_fasta = tempfile.NamedTemporaryFile()

    log.debug(f"Writing sequences to {input_fasta.name}")
    with open(input_fasta.name, "w") as f:
        for sequence in sequences:
            f.write(f">{sequence.id}\n{sequence.seq}\n")

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

    sequence_dict = {sequence.id: sequence for sequence in sequences}

    for record in SearchIO.parse(output_file.name, "hmmsearch3-domtab"):
        accession = record.accession.split('.')[0]
        record_name = accession_families[accession]

        for hit in record:
            for hsp in hit:
                # Skip RPW8 hits with low bitscores
                if record_name == "RPW8" and hsp.bitscore < 20:
                    continue

                # Lookup sequence by ID in the dictionary
                sequence = sequence_dict.get(hit.id)
                if sequence:
                    sequence.add_annotation(
                        record_name,
                        hsp.env_start,
                        hsp.env_end,
                        hsp.evalue,
                        hsp.bitscore,
                        "HMM",
                    )

    return sequences

