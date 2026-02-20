import pyhmmer
import pickle
import numpy as np
import os
from pathlib import Path
import logging
import warnings
from numpy.lib.stride_tricks import sliding_window_view
from resistify.annotation import Protein, Annotation

logger = logging.getLogger(__name__)

try:
    from sklearn.exceptions import InconsistentVersionWarning # type: ignore
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

MOTIF_MODELS = {
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

def load_models(search_type: str):
    models = {}
    base_path = Path(os.path.dirname(__file__)) / "data" / "nlrexpress_models"
    
    if search_type == "all":
        for predictor, path in MOTIF_MODELS.items():
            model_path = base_path / path
            models[predictor] = pickle.load(open(model_path, "rb"))
    elif search_type == "lrr":
        model_path = base_path / MOTIF_MODELS["LxxLxL"]
        models["LxxLxL"] = pickle.load(open(model_path, "rb"))
    return models

def nlrexpress(proteins: dict[str, Protein], search_type: str = "all", threads: int = 0):
    logger.info(f"Running NLRexpress to identify {search_type} motifs")
    models = load_models(search_type)
    db_path = Path(os.path.dirname(__file__)) / "data" / "nlrexpress.fa"

    alphabet = pyhmmer.easel.Alphabet.amino()
    queries = [
        pyhmmer.easel.DigitalSequence(
            name=p.id.encode(),
            sequence=alphabet.encode(p.sequence),
            alphabet=alphabet,
        ) for p in proteins.values()
    ]

    sequences = pyhmmer.easel.SequenceFile(db_path, digital=True, alphabet=alphabet)

    # Jackhmmer execution
    results = pyhmmer.hmmer.jackhmmer(
        queries,
        sequences,
        max_iterations=2,
        E=1e-5,
        domE=1e-5,
        checkpoints=True,
        cpus=threads,
    ) # type: ignore

    for result in results:
        sequence_id = result[0].hmm.name
        
        if len(result) < 2:
            logger.warning(f"No second iteration result for {sequence_id} - skipping")
            continue

        m1 = np.asarray(result[0].hmm.match_emissions[1:])
        m2 = np.asarray(result[1].hmm.match_emissions[1:])
        
        emission_matrix = -np.log(np.concatenate((m1, m2), axis=1))

        for motif_type, model in models.items():
            span = MOTIF_SPAN_LENGTHS[motif_type]
            window_size = 5 + span + 5 + 1
            
            if len(emission_matrix) < window_size:
                continue

            windows = sliding_window_view(emission_matrix, (window_size, 40))
            
            input_array = windows.reshape(-1, window_size * 40)

            predictions = model.predict_proba(input_array)
            
            hit_indices = np.where(predictions[:, 1] > 0.8)[0]

            for idx in hit_indices:
                score = predictions[idx, 1]
                proteins[sequence_id].add_annotation(
                    Annotation(
                        name=motif_type,
                        type="motif",
                        start=int(idx + 1),
                        end=int(idx + span),
                        source="nlrexpress",
                        score=float(score),
                    )
                )

    return proteins