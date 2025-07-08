import subprocess
import sys
import numpy as np
import pickle
from pathlib import Path
import tempfile
from multiprocessing import get_context
from concurrent.futures import ProcessPoolExecutor
from threadpoolctl import threadpool_limits
import shutil
import warnings
from resistify.utility import ProgressLogger, save_fasta
from resistify._loguru import logger
from resistify.annotations import Sequence

# Version 1.3 of sklearn introduced InconsistentVersionWarning, fall back to UserWarning if not available
# Necessary to suppress pickle version warnings
try:
    from sklearn.exceptions import InconsistentVersionWarning

    warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
except ImportError:
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

JACKHMMER_HEADER = [
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


def parse_jackhmmer(iteration_1_path: Path, iteration_2_path: Path):
    """
    Parse jackhmmer output for one or two iterations and return jackhmmer_results dict.
    If iteration_2_path is not provided or does not exist, only iteration 1 is used.
    """

    def parse_matrix(file: Path, iteration: bool = False):
        hmm_dict = {}
        with open(file) as f:
            i = 1
            lines = f.readlines()
            while i < len(lines):
                line = lines[i].split()
                if line[0] == "//":
                    i += 2
                elif line[0] == "NAME":
                    name = line[1]
                    if iteration:
                        name = name[:-3]
                    hmm_dict[name] = []
                    i += 15 if iteration else 14
                elif line[0] == "HMM":
                    i += 5
                else:
                    for value in line[1:21]:
                        try:
                            float(value)
                        except ValueError:
                            logger.error(f"{value} is not a float. {name} {file}")
                            sys.exit(1)
                    value_dict = dict(zip(JACKHMMER_HEADER, line[1:21]))
                    hmm_dict[name].append(value_dict)
                    i += 3
        return hmm_dict

    iteration_1 = parse_matrix(iteration_1_path, iteration=False)
    if iteration_2_path.exists():
        iteration_2 = parse_matrix(iteration_2_path, iteration=True)
    else:
        logger.debug(f"{iteration_2_path} does not exist, falling back to iteration 1")
        iteration_2 = iteration_1

    jackhmmer_results = {}
    for seq_id, matrix1 in iteration_1.items():
        jackhmmer_data = []
        for i, aa_dict in enumerate(matrix1):
            jackhmmer_data.append([])
            # Add iteration 1 values
            jackhmmer_data[-1].extend(list(aa_dict.values()))
            # Add iteration 2 values if available, else repeat iteration 1
            if iteration_2 and seq_id in iteration_2:
                jackhmmer_data[-1].extend(list(iteration_2[seq_id][i].values()))
            else:
                jackhmmer_data[-1].extend(list(aa_dict.values()))
        jackhmmer_results[seq_id] = jackhmmer_data
    return jackhmmer_results


class MotifPredictor:
    def __init__(self):
        self.motif_models = {
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
        self.MOTIF_SPAN_LENGTHS = {
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
        self.models = {}

    def load(self, search_type):
        self.models = {}
        if search_type == "all":
            for predictor, path in self.motif_models.items():
                model_path = Path(__file__).parent / "data" / path
                logger.debug(f"Loading model for {predictor}")
                model = pickle.load(open(model_path, "rb"))
                self.models[predictor] = model
        elif search_type == "lrr":
            model_path = Path(__file__).parent / "data" / self.motif_models["LxxLxL"]
            model = pickle.load(open(model_path, "rb"))
            self.models["LxxLxL"] = model
        return self.models

    def predict(self, sequences, jackhmmer_results):
        for predictor, model in self.models.items():
            logger.debug(f"Predicting {predictor} for batch")
            motif_size = self.MOTIF_SPAN_LENGTHS[predictor]
            matrix = []
            for sequence in sequences:
                sequence_length = len(sequence.seq)
                for i in range(sequence_length):
                    if i >= 5 and i < sequence_length - (motif_size + 5):
                        matrix.append([])
                        for j in range(-5, motif_size + 6):
                            matrix[-1] += jackhmmer_results[sequence.id][i + j]
            if not matrix:
                continue
            matrix = np.array(matrix, dtype=float)
            with threadpool_limits(limits=2):
                result = model.predict_proba(matrix)
            result_index = 0
            for sequence in sequences:
                sequence_length = len(sequence.seq)
                for i in range(sequence_length):
                    if i >= 5 and i < sequence_length - (motif_size + 5):
                        probability = round(result[result_index][1], 4)
                        if probability > 0.8:
                            sequence.add_motif(predictor, probability, i, motif_size)
                        result_index += 1
        return sequences


def nlrexpress(
    sequences: list[Sequence], search_type, chunk_size, threads: int
) -> list[Sequence]:
    # Test that jackhmmer is functional
    try:
        subprocess.run(
            ["jackhmmer", "-h"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    except Exception as e:
        logger.critical(e.stderr.strip() if e.stderr else e)
        sys.exit(1)

    predictor = MotifPredictor()
    predictor.load(search_type)

    jackhmmer_db = tempfile.NamedTemporaryFile()
    nlrexpress_database = Path(__file__).parent / "data" / "nlrexpress.fasta"
    shutil.copyfile(nlrexpress_database, jackhmmer_db.name)

    batches = [
        sequences[i : i + chunk_size] for i in range(0, len(sequences), chunk_size)
    ]

    logger.info("Running NLRexpress - this could take a while...")

    progress_logger = ProgressLogger(len(batches))

    executor = ProcessPoolExecutor(
        max_workers=-(-threads // 2), mp_context=get_context("spawn")
    )

    process_arguments = [(batch, jackhmmer_db.name, predictor) for batch in batches]

    results = []

    try:
        for result in executor.map(nlrexpress_subprocess, process_arguments):
            progress_logger.update()
            results.extend(result)
    except KeyboardInterrupt:
        executor.shutdown(wait=False, cancel_futures=True)
        sys.exit(1)

    if len(results) != len(sequences):
        logger.critical(
            "Sequences dropped during NLRexpress - this should not happen and must be reported"
        )
        sys.exit(1)

    return results


@logger.catch
def nlrexpress_subprocess(params):
    sequences, jackhmmer_db, predictor = params

    temp_dir = tempfile.TemporaryDirectory(delete=False)
    fasta_path = Path(temp_dir.name) / "seq.fa"

    save_fasta(sequences, fasta_path)

    out_path = fasta_path.parent / (fasta_path.name + ".out")
    iteration_1_path = fasta_path.parent / (fasta_path.name + ".out-1.hmm")
    iteration_2_path = fasta_path.parent / (fasta_path.name + ".out-2.hmm")

    cmd = [
        "jackhmmer",
        "--noali",
        "-N",
        "2",
        "--cpu",
        "2",
        "-E",
        "1e-5",
        "--domE",
        "1e-5",
        "--chkhmm",
        str(out_path),
        str(fasta_path),
        str(jackhmmer_db),
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
        logger.critical(f"stdout: {e.stderr.strip()}")
        logger.critical(f"stderr: {e.stderr.strip()}")
        raise

    logger.debug(f"Jackhmmer for {fasta_path} has completed")

    jackhmmer_results = parse_jackhmmer(iteration_1_path, iteration_2_path)

    predictor.predict(sequences, jackhmmer_results)

    temp_dir.cleanup()

    logger.debug(f"NLRexpress subprocess for {fasta_path} has completed")

    return sequences
