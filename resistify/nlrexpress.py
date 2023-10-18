import subprocess
import sys
import numpy as np
from sklearn.neural_network import MLPClassifier
import pickle
from Bio import SeqIO
import tempfile
import logging
import os
from resistify.annotations import Motif

motif_span_lengths = {
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
}


def jackhmmer(input_fasta, sequences, temp_dir, data_dir):
    """
    Run jackhmmer on the input fasta file against the database_path.
    """

    database_file = os.path.join(data_dir, "nlrexpress.fasta")

    cmd = [
        "jackhmmer",
        "--noali",
        "-N",
        "2",  # number of iterations
        "--cpu",
        "4",
        "-E",
        "1e-5",
        "--domE",
        "1e-5",
        "--chkhmm",
        temp_dir.name + "/jackhmmer",
        input_fasta,
        database_file,
    ]

    logging.info(f"😊 Running jackhmmer...")
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        logging.info(f"😊 jackhmmer completed successfully...")
    except subprocess.CalledProcessError as e:
        logging.error(f"😞 Error running jackhmmer. Stdout of jackhmmer: {e.stdout}")
        sys.exit(1)

    jackhmmer_iteration_1 = parse_jackhmmer(
        os.path.join(temp_dir.name, "jackhmmer-1.hmm"), iteration=False
    )
    jackhmmer_iteration_2 = parse_jackhmmer(
        os.path.join(temp_dir.name, "jackhmmer-2.hmm"), iteration=True
    )

    sequences = prepare_jackhmmer_data(
        sequences, jackhmmer_iteration_1, jackhmmer_iteration_2
    )

    return sequences


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
                        print(f"Error: {value} is not a float. {name} {file}")
                        sys.exit(1)
                value_dict = dict(zip(header1, line[1:21]))
                hmm_dict[name].append(value_dict)
                i += 3

    return hmm_dict


def parse_fasta(file):
    """
    Parse a FASTA file and return a dictionary of the sequences.
    """
    seq_dict = {}

    logging.info(f"😊 Parsing FASTA file...")
    with open(file) as f:
        for record in SeqIO.parse(f, "fasta"):
            seq_dict[record.id] = record.seq

    return seq_dict


def prepare_jackhmmer_data(sequences, hmm_it1, hmm_it2):
    for sequence in sequences:
        # get the dna sequence
        seq = sequences[sequence].sequence
        # make blank list for this sequence
        jackhmmer_data = []

        # for each amino acid in the sequence
        for i, aa in enumerate(seq):
            # add a blank list for this amino acid
            jackhmmer_data.append([])
            for k, (key, val) in enumerate(hmm_it1[sequence][i].items()):
                jackhmmer_data[-1].append(val)

            if sequence in hmm_it2:
                for k, (key, val) in enumerate(hmm_it2[sequence][i].items()):
                    jackhmmer_data[-1].append(val)
            else:
                for k, (key, val) in enumerate(hmm_it1[sequence][i].items()):
                    jackhmmer_data[-1].append(val)

        sequences[sequence].jackhmmer_data = jackhmmer_data

    return sequences


def predict_motif(sequences, predictor, data_dir):
    # create a matrix for the predictors motif size
    logging.info(f"😊 Generating matrix for {predictor}...")
    matrix = []
    motif_size = motif_span_lengths[predictor]
    for sequence in sequences:
        sequence_length = len(sequences[sequence].sequence)
        for i in range(sequence_length):
            # make sure we are within the sequence bounds
            if i >= 5 and i < sequence_length - (motif_size + 5):
                features = sequences[sequence].jackhmmer_data

                matrix.append([])

                # add the hmm values for this range
                for j in range(-5, motif_size + 6):
                    matrix[-1] += features[i + j]

    matrix = np.array(matrix, dtype=float)
    # load the model from model dictionary

    model_path = os.path.join(data_dir, motif_models[predictor])

    model = pickle.load(open(model_path, "rb"))
    # run the prediction
    logging.info(f"😊 Predicting {predictor} motifs...")
    result = model.predict_proba(matrix)

    # iterate through sequences and pull out any predictions with a probability > 0.8
    result_index = 0
    for sequence in sequences:
        sequence_length = len(sequences[sequence].sequence)
        for i in range(len(sequences[sequence].sequence)):
            # make sure we are within the sequence bounds
            if i >= 5 and i < sequence_length - (motif_span_lengths[predictor] + 5):
                value = round(result[result_index][1], 4)
                if value > 0.8:
                    sequences[sequence].add_motif(Motif(predictor, value, i))
                result_index += 1
