import sys
import warnings
import numpy as np
import pickle
import os
import logging
import pyhmmer
import threadpoolctl
from tqdm.auto import tqdm
from resistify.annotation import Annotation, Protein

# Version 1.3 of sklearn introduced InconsistentVersionWarning, fall back to UserWarning if not available
# Necessary to suppress pickle version warnings
try:
    from sklearn.exceptions import InconsistentVersionWarning  # type: ignore

    warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
except ImportError:
    warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

logger = logging.getLogger(__name__)

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
    "RNBS-A": 10,
    "Walker-B": 8,
    "RNBS-B": 7,
    "RNBS-C": 10,
    "RNBS-D": 9,
    "GLPL": 5,
    "MHD": 3,
    "LxxLxL": 6,
}

MOTIF_MODELS = {
    "extEDVID": "MLP_CC_extEDVID.pkl",
    "VG": "MLP_NBS_VG.pkl",
    "P-loop": "MLP_NBS_P-loop.pkl",
    "RNBS-A": "MLP_NBS_RNBS-A.pkl",
    "RNBS-B": "MLP_NBS_RNBS-B.pkl",
    "RNBS-C": "MLP_NBS_RNBS-C.pkl",
    "RNBS-D": "MLP_NBS_RNBS-D.pkl",
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


def _load_models(search_type: str):
    candidates = (
        MOTIF_MODELS
        if search_type == "all"
        else {search_type: MOTIF_MODELS[search_type]}
    )
    models = {}
    for predictor, filename in candidates.items():
        path = os.path.join(
            os.path.dirname(__file__), "data", "nlrexpress_models", filename
        )
        with open(path, "rb") as f:
            models[predictor] = pickle.load(f)
    return models


def _predict_motifs(
    protein: Protein, mat1: np.ndarray, mat2: np.ndarray, models: dict[str, object]
):
    seq_len = protein.length
    mat_combined = np.concatenate([mat1[1:], mat2[1:]], axis=1)

    for predictor, model in models.items():
        motif_size = MOTIF_SPAN_LENGTHS[predictor]
        window_size = motif_size + 11

        if seq_len < window_size:
            continue

        windows = np.lib.stride_tricks.sliding_window_view(
            mat_combined, window_size, axis=0
        )
        matrix = windows.transpose(0, 2, 1).reshape(windows.shape[0], -1)

        proba = model.predict_proba(matrix)

        for k, prob in enumerate(proba):
            if prob[1] > 0.8:
                i = k + 5
                protein.add_annotation(
                    Annotation(
                        name=predictor,
                        start=i + 1,
                        end=i + motif_size,
                        type="motif",
                        source="nlrexpress",
                        score=round(float(prob[1]), 4),
                    )
                )


def nlrexpress(proteins: dict[str, Protein], threads: int, search_type: str = "all"):

    models = _load_models(search_type)

    alphabet = pyhmmer.easel.Alphabet.amino()
    db_path = os.path.join(os.path.dirname(__file__), "data", "nlrexpress.fasta")

    queries = [
        pyhmmer.easel.TextSequence(
            name=protein.id.encode(), sequence=protein.sequence
        ).digitize(alphabet)
        for protein in proteins.values()
    ]

    with pyhmmer.easel.SequenceFile(db_path, digital=True, alphabet=alphabet) as f:
        database = f.read_block()

    logger.info("Running NLRexpress...")

    progress = tqdm(
        total=len(queries), desc="NLRexpress", disable=not sys.stdout.isatty()
    )

    def _progress_callback(query, total):
        progress.update(1)

    with threadpoolctl.threadpool_limits(limits=threads):
        for iteration_results in pyhmmer.hmmer.jackhmmer(
            queries,
            database,
            max_iterations=2,
            checkpoints=True,
            cpus=threads,
            callback=_progress_callback,
            E=1e-5,
            domE=1e-5,
        ):
            seq_id = iteration_results[0].hmm.name
            if isinstance(seq_id, bytes):
                seq_id = seq_id.decode()

            mat1 = -np.log(np.array(iteration_results[0].hmm.match_emissions) + 1e-8)

            if len(iteration_results) == 2:
                mat2 = -np.log(
                    np.array(iteration_results[1].hmm.match_emissions) + 1e-8
                )
            else:
                mat2 = mat1

            _predict_motifs(proteins[seq_id], mat1, mat2, models)

    progress.close()

    return proteins
