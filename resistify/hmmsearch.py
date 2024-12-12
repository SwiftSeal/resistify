import subprocess
import sys
import os
import logging
import tempfile

log = logging.getLogger(__name__)

accession_families = {
    "PF00931": ("NB-ARC", 23.5),
    "PF01582": ("TIR", 21.3),
    "PF05659": ("RPW8", 30.4),
    "PF13676": ("TIR", 28),
    "PF18052": ("CC", 27.7),
    "PF20160": ("C-JID", 24.1),
    "MADA": ("MADA", 10),  # From original paper
    "PF00069": ("PKinase", 31.7),
    "PF01453": ("G-LecRLK", 32.4),
    "PF00954": ("G-LecRLK", 24),
    "PF08276": ("G-LecRLK", 21.1),
    "PF00024": ("G-LecRLK", 21.4),
    "PF14295": ("G-LecRLK", 27.1),
    "PF00139": ("L-LecRLK", 26.3),
    "PF00059": ("C-LecRLK", 22.6),
    "PF13947": ("WAK", 29.4),
    "PF14380": ("WAK", 27),
    "PF00008": ("WAK", 21.5),
    "PF08488": ("WAK", 22.1),
    "PF07645": ("WAK", 27),
    "PF12662": ("WAK", 26.6),
    "PF12947": ("WAK", 28.5),
    "PF11721": ("CrRLK1L", 23),
    "PF12819": ("CrRLK1L", 24.5),
    "PF01476": ("LysM", 20.9),
    "PF01657": ("CRK", 25),
    "PF00314": ("Thaumatin", 27.4),
    "PF13540": ("CR-like", 22),
    "PF19160": ("SPARK", 25),
    "PF00704": ("GH18", 29.6),
    "PF00182": ("GH19", 25.5),
    "PF00188": ("CAP", 21.1),
    "PF16101": ("PRIMA1", 38.6),
}


def hmmsearch(sequences, search_type):
    if search_type == "nlr":
        hmmsearch_db = os.path.join(os.path.dirname(__file__), "data", "nlrdb.hmm")
    elif search_type == "prr":
        hmmsearch_db = os.path.join(os.path.dirname(__file__), "data", "prrdb.hmm")
    output_file = tempfile.NamedTemporaryFile()

    input_fasta = tempfile.NamedTemporaryFile()

    log.debug(f"Writing sequences to {input_fasta.name}")
    with open(input_fasta.name, "w") as f:
        for sequence in sequences:
            f.write(f">{sequence.id}\n{sequence.seq}\n")

    cmd = [
        "hmmsearch",
        "--noali",
        "-Z 45638612",
        # "--domE",
        # evalue,
        "--domtblout",
        output_file.name,
        hmmsearch_db,
        input_fasta.name,
    ]

    try:
        log.info("Running hmmsearch...")
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        log.error(f"Error running hmmsearch:\nStderr: {e.stderr}\nStdout:{e.stdout}")
        sys.exit(1)
    except FileNotFoundError:
        log.error("hmmsearch not found. Have you installed it?")
        sys.exit(1)

    sequence_dict = {sequence.id: sequence for sequence in sequences}

    with open(output_file.name) as file:
        for line in file:
            if line.startswith("#"):
                continue
            line = line.split()
            sequence_id = line[0]
            accession = line[4].split(".")[0]
            evalue = float(line[11])
            bitscore = float(line[13])
            start = int(line[17])
            end = int(line[18])

            domain_name, domain_bit_threshold = accession_families[accession]

            if bitscore < domain_bit_threshold:
                    continue

                # Lookup sequence by ID in the dictionary
            sequence = sequence_dict.get(sequence_id)
            if sequence:
                sequence.add_annotation(
                    domain_name,
                    "HMM",
                    start,
                    end,
                    evalue,
                    bitscore,
                )

    return sequences
