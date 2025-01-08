import subprocess
import sys
import numpy as np
import pickle
import os
import tempfile
from multiprocessing import Pool, cpu_count, get_context
from threadpoolctl import threadpool_limits
import shutil
import warnings
from resistify.utility import log_percentage
from resistify._loguru import logger

# Version 1.3 of sklearn introduced InconsistentVersionWarning, fall back to UserWarning if not available
# Necessary to suppress pickle version warnings
try:
    from sklearn.exceptions import InconsistentVersionWarning

    warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
except ImportError:
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

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


def parse_jackhmmer(file, iteration=False):
    logger.debug(f"Parsing jackhmmer {file} as iteration {iteration}")
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
                if not iteration:
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
                        logger.error(f"{value} is not a float. {name} {file}")
                        sys.exit(1)
                value_dict = dict(zip(header1, line[1:21]))
                hmm_dict[name].append(value_dict)
                i += 3

    return hmm_dict


def nlrexpress(sequences, search_type, chunk_size, threads):
    total_sequences = len(sequences)
    if threads is None:
        try:
            threads = len(os.sched_getaffinity(0))
        except AttributeError:
            threads = cpu_count()

    models = load_models(search_type)

    with tempfile.NamedTemporaryFile(delete=False) as jackhmmer_db:
        nlrexpress_database = os.path.join(
            os.path.dirname(__file__), "data", "nlrexpress.fasta"
        )
        shutil.copyfile(nlrexpress_database, jackhmmer_db.name)

    batches = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]

    args = [(batch, jackhmmer_db.name, models) for batch in batches]

    results = []
    logger.info("Running NLRexpress - this could take a while...")

    total_iterations = len(args)
    iteration = 0
    # Need to use spawn otherwise the subprocesses will hang
    with get_context("spawn").Pool(-(-threads // 2)) as pool:
        for result in pool.imap(nlrexpress_subprocess, args):
            results.append(result)
            iteration += 1
            log_percentage(iteration, total_iterations)

    sequences = [seq for batch in results for seq in batch]

    if len(sequences) != total_sequences:
        log.critical("Sequences dropped during NLRexpress - this should not happen and must be reported")

    return sequences


def load_models(search_type):
    models = {}
    if search_type == "all":
        for predictor, path in motif_models.items():
            model_path = os.path.join(os.path.dirname(__file__), "data", path)
            logger.debug(f"Loading model for {predictor}")
            model = pickle.load(open(model_path, "rb"))
            models[predictor] = model
    elif search_type == "lrr":
        model_path = os.path.join(
            os.path.dirname(__file__), "data", motif_models["LxxLxL"]
        )
        model = pickle.load(open(model_path, "rb"))
        models["LxxLxL"] = model
    return models


def nlrexpress_subprocess(params):
    sequences, jackhmmer_db, models = params

    temp_dir = tempfile.TemporaryDirectory()
    fasta_path = os.path.join(temp_dir.name, "seq.fa")
    iteration_1_path = fasta_path + ".out-1.hmm"
    iteration_2_path = fasta_path + ".out-2.hmm"

    with open(fasta_path, "w") as f:
        for sequence in sequences:
            logger.debug(f"Writing {sequence.id} to {fasta_path}")
            f.write(f">{sequence.id}\n{sequence.seq}\n")

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
        fasta_path + ".out",
        fasta_path,
        jackhmmer_db,
    ]

    try:
        logger.debug(f"Running jackhmmer for {fasta_path}")
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except subprocess.CalledProcessError as e:
        sys.exit(1)

    iteration_1 = parse_jackhmmer(iteration_1_path, iteration=False)

    if os.path.exists(iteration_2_path):
        iteration_2 = parse_jackhmmer(iteration_2_path, iteration=True)
    else:
        iteration_2 = parse_jackhmmer(iteration_1_path, iteration=False)

    jackhmmer_results = {}

    for sequence in sequences:
        # make blank list for this sequence
        jackhmmer_data = []

        # for each amino acid in the sequence
        for i, aa in enumerate(sequence.seq):
            # add a blank list for this amino acid
            jackhmmer_data.append([])
            for k, (key, val) in enumerate(iteration_1[sequence.id][i].items()):
                jackhmmer_data[-1].append(val)

            if sequence.id in iteration_2:
                for k, (key, val) in enumerate(iteration_2[sequence.id][i].items()):
                    jackhmmer_data[-1].append(val)
            else:
                for k, (key, val) in enumerate(iteration_1[sequence.id][i].items()):
                    jackhmmer_data[-1].append(val)

        jackhmmer_results[sequence.id] = jackhmmer_data

    for predictor, model in models.items():
        logger.debug(f"Predicting {predictor} for result of {fasta_path}")
        motif_size = MOTIF_SPAN_LENGTHS[predictor]
        matrix = []
        for sequence in sequences:
            sequence_length = len(sequence.seq)
            for i in range(sequence_length):
                if i >= 5 and i < sequence_length - (motif_size + 5):
                    matrix.append([])
                    for j in range(-5, motif_size + 6):
                        matrix[-1] += jackhmmer_results[sequence.id][i + j]

        matrix = np.array(matrix, dtype=float)

        with threadpool_limits(limits=2):
            result = model.predict_proba(matrix)

        result_index = 0
        for sequence in sequences:
            sequence_length = len(sequence.seq)
            for i in range(sequence_length):
                if i >= 5 and i < sequence_length - (motif_size + 5):
                    value = round(result[result_index][1], 4)
                    if value > 0.8:
                        sequence.add_motif(predictor, value, i)
                    result_index += 1

    temp_dir.cleanup()

    logger.debug(f"NLRexpress subprocess for {fasta_path} has completed")

    return sequences
