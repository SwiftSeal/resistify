import subprocess
import sys
import numpy as np
from sklearn.neural_network import MLPClassifier
import pickle
from Bio import SeqIO
import tempfile
import logging

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
        'VG': 'models/MLP_NBS_VG.pkl',
        'P-loop': 'models/MLP_NBS_P-loop.pkl',
        'RNSB-A': 'models/MLP_NBS_RNSB-A.pkl',
        'RNSB-B': 'models/MLP_NBS_RNSB-B.pkl',
        'RNSB-C': 'models/MLP_NBS_RNSB-C.pkl',
        'RNSB-D': 'models/MLP_NBS_RNSB-D.pkl',
        'Walker-B': 'models/MLP_NBS_Walker-B.pkl',
        'GLPL': 'models/MLP_NBS_GLPL.pkl',
        'MHD': 'models/MLP_NBS_MHD.pkl'
    }

def jackhmmer(input_fasta, temp_dir, database_path):
    """
    Run jackhmmer on the input fasta file against the database_path.
    """

    cmd = [
        "jackhmmer",
        "--cpu",
        "6",
        "-N",
        "2", # number of iterations
        "-E",
        "1e-5",
        "--domE",
        "1e-5",
        "--chkhmm",
        temp_dir.name + "/jackhmmer",
        input_fasta,
        database_path,
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
    

def parse_jackhmmer(file, iteration = False):
    hmm_dict = {}

    header1 = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

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

def generateInputFile(sequences, hmm_it1, hmm_it2):
    data = {}

    for name in sequences:
        # get the dna sequence
        seq = sequences[name]
        # make blank list for this sequence
        data[name] = []

        # for each amino acid in the sequence
        for i, aa in enumerate(seq):
            # add a blank list for this amino acid
            data[name].append([])
            for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                data[name][-1].append(val)

            if name in hmm_it2:
                for k, (key, val) in enumerate( hmm_it2[name][i].items() ):
                    data[name][-1].append(val)
            else:
                for k, (key, val) in enumerate( hmm_it1[name][i].items() ):
                    data[name][-1].append(val)

    return data

def generate_matrix(sequences, input_data, predictor):
    logging.info(f"😊 Generating matrix for {predictor}...")
    matrix = []
    motif_size = motif_span_lengths[predictor]
    for sequence in sequences:
        sequence_length = len(sequences[sequence])
        for i in range(sequence_length):
            # make sure we are within the sequence bounds
            if i >= 5 and i < sequence_length - (motif_size + 5):
                features = input_data[sequence]

                matrix.append([])

                # add the hmm values for this range
                for j in range(-5, motif_size + 6):
                    matrix[-1] += features[i + j]

    matrix = np.array(matrix, dtype=float)
    return matrix

def predict_motif(sequences, input_data, predictor):
    # create a matrix for the predictors motif size
    matrix = generate_matrix(sequences, input_data, predictor)
    # load the model from model dictionary
    model = pickle.load(open(motif_models[predictor], 'rb'))
    # run the prediction
    result = model.predict_proba(matrix)

    return result
