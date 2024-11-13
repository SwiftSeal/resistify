import subprocess
import sys
import numpy as np
import pickle
import os
import logging
import tempfile
from sklearn.neural_network import MLPClassifier
from multiprocessing import Pool
from resistify.annotations import Motif
import shutil

log = logging.getLogger(__name__)

MOTIF_SPAN_LENGTHS = {
    "extEDVID": 12,
    "bA": 10,
    "aA": 7,
    "bC": 8,
    "aC": 6,
    "bDaD1": 16,
    "aD3": 13,
    "VG": 5,
    "P-loop": 9,
    "RNSB-A": 10,
    "Walker-B": 8,
    "RNSB-B": 7,
    "RNSB-C": 10,
    "RNSB-D": 9,
    "GLPL": 5,
    "MHD": 3,
    "LxxLxL": 6,
}

motif_models = {
    "extEDVID": "MLP_CC_extEDVID.pkl",
    "VG": "MLP_NBS_VG.pkl",
    "P-loop": "MLP_NBS_P-loop.pkl",
    "RNSB-A": "MLP_NBS_RNSB-A.pkl",
    "RNSB-B": "MLP_NBS_RNSB-B.pkl",
    "RNSB-C": "MLP_NBS_RNSB-C.pkl",
    "RNSB-D": "MLP_NBS_RNSB-D.pkl",
    "Walker-B": "MLP_NBS_Walker-B.pkl",
    "GLPL": "MLP_NBS_GLPL.pkl",
    "MHD": "MLP_NBS_MHD.pkl",
    "LxxLxL": "MLP_LRR_LxxLxL.pkl",
    "aA": "MLP_TIR_aA.pkl",
    "aC": "MLP_TIR_aC.pkl",
    "aD3": "MLP_TIR_aD3.pkl",
    "bA": "MLP_TIR_bA.pkl",
    "bC": "MLP_TIR_bC.pkl",
    "bDaD1": "MLP_TIR_bD-aD1.pkl",
}


def split_fasta(sequences, chunk_size):
    """
    Split a fasta file into chunks of defined size.
    Return a list of the file paths.
    """
    log.debug(f"Splitting fasta into chunks of {chunk_size} sequences")
    fastas = []

    for i in range(0, len(sequences), chunk_size):
        chunk = sequences[i : i + chunk_size]
        chunk_path = tempfile.NamedTemporaryFile(delete=False)
        log.debug(f"Writing chunk to {chunk_path.name}")
        with open(chunk_path.name, "w") as f:
            for sequence in chunk:
                f.write(f">{sequence.id}\n{sequence.seq}\n")
        fastas.append(chunk_path.name)

    return fastas


def jackhmmer_subprocess(fasta, database):
    """
    Run jackhmmer on the input fasta file against the database_path.
    """
    cmd = [
        "jackhmmer",
        "--noali",
        "-N",
        "2",  # number of iterations
        "--cpu",
        "2",
        "-E",
        "1e-5",
        "--domE",
        "1e-5",
        "--chkhmm",
        fasta + ".out",
        fasta,
        database,
    ]
    try:
        log.debug(f"Running jackhmmer on {fasta}")
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        log.error(f"Error running jackhmmer:\nStderr: {e.stderr}\nStdout:{e.stdout}")
        sys.exit(1)


def jackhmmer(fastas, threads):
    """
    Run jackhmmer on the input fasta file against the database_path.
    """
    # Moving jackhmmer to temp can help speed stuff up on clusters with local storage.
    with tempfile.NamedTemporaryFile(delete=False) as temp_database:
        nlrexpress_database = os.path.join(
            os.path.dirname(__file__), "data", "nlrexpress.fasta"
        )
        shutil.copyfile(nlrexpress_database, temp_database.name)

        # run jackhmmer on each chunk
        log.info(f"Running jackhmmer, this could take a while...")
        args = [(fasta, temp_database.name) for fasta in fastas]
        with Pool(-(-threads // 2)) as pool:
            pool.starmap(jackhmmer_subprocess, args)

    # merge the chunks
    iteration_1_path = tempfile.NamedTemporaryFile()
    with open(iteration_1_path.name, "w") as f:
        for fasta in fastas:
            log.debug(f"Writing {fasta}.out-1.hmm to {iteration_1_path.name}")
            with open(f"{fasta}.out-1.hmm") as chunk:
                f.write(chunk.read())

    jackhmmer_iteration_1 = parse_jackhmmer(iteration_1_path.name, iteration=False)

    iteration_2_path = tempfile.NamedTemporaryFile()
    with open(iteration_2_path.name, "w") as f:
        for fasta in fastas:
            log.debug(f"Writing {fasta}.out-2.hmm to {iteration_2_path.name}")
            try:
                with open(f"{fasta}.out-2.hmm") as chunk:
                    f.write(chunk.read())
            except FileNotFoundError:
                log.debug("Second jackhmmer iteration file does not exist, skipping...")

    try:
        jackhmmer_iteration_2 = parse_jackhmmer(iteration_2_path.name, iteration=True)
    except FileNotFoundError:
        log.debug(
            f"Second jackhmmer iteration file does not exist, setting second as first..."
        ),
        jackhmmer_iteration_2 = parse_jackhmmer(iteration_1_path.name, iteration=False)

    return jackhmmer_iteration_1, jackhmmer_iteration_2


def parse_jackhmmer(file, iteration=False):
    hmm_dict = {}

    header1 = [
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    ]

    with open(file) as f:
        i = 1
        lines = f.readlines()
        while i < len(lines):
            line = lines[i].split()
            if line[0] == "//":
                # if we hit // then we are onto the new entry - skip forward 2 lines to get to next NAME
                i += 2
            elif line[0] == "NAME":
                name = line[1]
                # if this is an iteration, need to strip "-i<iteration number>" from the name
                if iteration:
                    name = name[:-3]
                hmm_dict[name] = []
                # There is an extra field in the iteration file, so we need to skip an extra line
                if iteration == False:
                    i += 14
                else:
                    i += 15
            elif line[0] == "HMM":
                # we have hit the start of the matrix, skip forward 5 lines to get to the matrix
                i += 5
            else:
                # we are in a matrix entry for an amino acid
                # create a dictionary of amino acids to values, then skip forward 3 lines to move to next amino acid
                # ensure that all values are floats
                for value in line[1:21]:
                    try:
                        float(value)
                    except ValueError:
                        log.error(f"{value} is not a float. {name} {file}")
                        sys.exit(1)
                value_dict = dict(zip(header1, line[1:21]))
                hmm_dict[name].append(value_dict)
                i += 3

    return hmm_dict


def prepare_jackhmmer_data(sequences, hmm_it1, hmm_it2):
    for sequence in sequences:
        log.debug(f"Preparing jackhmmer data for {sequence.id}")
        # make blank list for this sequence
        jackhmmer_data = []

        # for each amino acid in the sequence
        for i, aa in enumerate(sequence.seq):
            # add a blank list for this amino acid
            jackhmmer_data.append([])
            for k, (key, val) in enumerate(hmm_it1[sequence.id][i].items()):
                jackhmmer_data[-1].append(val)

            if sequence.id in hmm_it2:
                for k, (key, val) in enumerate(hmm_it2[sequence.id][i].items()):
                    jackhmmer_data[-1].append(val)
            else:
                for k, (key, val) in enumerate(hmm_it1[sequence.id][i].items()):
                    jackhmmer_data[-1].append(val)

        sequence.jackhmmer_data = jackhmmer_data

    return sequences


def predict_motif(sequences, predictor):
    log.info(f"Predicting {predictor} motifs...")
    motif_size = MOTIF_SPAN_LENGTHS[predictor]
    total_length = sum(len(sequence.seq) - (motif_size + 10) for sequence in sequences)
    
    # Preallocate the matrix with the correct size
    matrix = np.zeros((total_length, (motif_size + 11) * len(sequences[0].jackhmmer_data[0])), dtype=float)
    
    index = 0
    for sequence in sequences:
        log.debug(f"Preparing matrix for {sequence.id}")
        features = sequence.jackhmmer_data
        for i in range(5, len(sequence.seq) - (motif_size + 5)):
            matrix[index] = np.concatenate([features[i + j] for j in range(-5, motif_size + 6)])
            index += 1

    # Load the model from model dictionary
    model_path = os.path.join(os.path.dirname(__file__), "data", motif_models[predictor])
    model = pickle.load(open(model_path, "rb"))
    log.debug("Running prediction")
    result = model.predict_proba(matrix)
    log.debug("Prediction complete")

    result_index = 0
    for sequence in sequences:
        log.debug(f"Adding motifs to {sequence.id}")
        for i in range(5, len(sequence.seq) - (motif_size + 5)):
            value = round(result[result_index][1], 4)
            if value > 0.8:
                sequence.add_motif(Motif(predictor, value, i))
            result_index += 1

    log.info(f"Finished predicting {predictor} motifs.")


def nlrexpress(sequences, chunk_size, threads):

    fastas = split_fasta(sequences, chunk_size)

    jackhmmer_iteration_1, jackhmmer_iteration_2 = jackhmmer(fastas, threads)

    sequences = prepare_jackhmmer_data(
        sequences, jackhmmer_iteration_1, jackhmmer_iteration_2
    )

    for predictor in motif_models.keys():
        predict_motif(sequences, predictor)

    return sequences
