import pyhmmer
import pickle
import numpy as np
import os
from pathlib import Path
import logging
import warnings
from tqdm.auto import tqdm
from resistify.annotation import Protein, Annotation

logger = logging.getLogger(__name__)

# Version 1.3 of sklearn introduced InconsistentVersionWarning, fall back to UserWarning if not available
# Necessary to suppress pickle version warnings
# This will probably blow up at some point in the future
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
    if search_type == "all":
        for predictor, path in MOTIF_MODELS.items():
            model_path = (
                Path(os.path.dirname(__file__)) / "data" / "nlrexpress_models" / path
            )
            logger.debug(f"Loading model for {predictor}")
            model = pickle.load(open(model_path, "rb"))
            models[predictor] = model
    elif search_type == "lrr":
        model_path = (
            Path(os.path.dirname(__file__))
            / "data"
            / "nlrexpress_models"
            / MOTIF_MODELS["LxxLxL"]
        )
        model = pickle.load(open(model_path, "rb"))
        models["LxxLxL"] = model
    return models


def nlrexpress(
    proteins: dict[str, Protein], search_type: str = "all", threads: int = 0
):
    logger.info(f"Running NLRexpress to identify {search_type} motifs")
    models = load_models(search_type)
    db_path = Path(os.path.dirname(__file__)) / "data" / "nlrexpress.fa"

    alphabet = pyhmmer.easel.Alphabet.amino()

    queries = []
    for protein in proteins.values():
        queries.append(
            pyhmmer.easel.DigitalSequence(
                name=protein.id.encode(),
                sequence=alphabet.encode(protein.sequence),
                alphabet=alphabet,
            )
        )

    sequences = pyhmmer.easel.SequenceFile(db_path, digital=True, alphabet=alphabet)

    with tqdm(total=len(proteins)) as pbar:
        for result in pyhmmer.hmmer.jackhmmer(
            queries,
            sequences,
            max_iterations=2,
            E=1e-5,
            domE=1e-5,
            checkpoints=True,
            cpus=threads,
            callback=lambda _, __: pbar.update(1),
        ):
            sequence_id = result[0].hmm.name
            logger.debug(f"Processing {sequence_id}")
            try:
                emission_matrix = np.concatenate(
                    (
                        np.asarray(result[0].hmm.match_emissions[1:]),
                        np.asarray(result[1].hmm.match_emissions[1:]),
                    ),
                    axis=1,
                )
            except IndexError:
                # This will occur if the second iteration fails.
                logger.warning(
                    f"No second iteration result for {sequence_id} - skipping"
                )
                continue

            emission_matrix = -np.log(emission_matrix)

            for motif_type, model in models.items():
                logger.debug(f"Searching for {motif_type} in {sequence_id}")
                input_rows = []
                for i in range(
                    5, len(emission_matrix) - MOTIF_SPAN_LENGTHS[motif_type] - 5
                ):
                    flattened = emission_matrix[
                        i - 5 : i + MOTIF_SPAN_LENGTHS[motif_type] + 5 + 1
                    ].flatten()
                    input_rows.append(flattened)

                input_array = np.array(input_rows)

                predictions = model.predict_proba(input_array)

                for i, prob in enumerate(predictions):
                    if prob[1] > 0.8:
                        proteins[sequence_id].add_annotation(
                            Annotation(
                                name=motif_type,
                                type="motif",
                                start=i + 1,
                                end=i + MOTIF_SPAN_LENGTHS[motif_type],
                                source="nlrexpress",
                                score=prob[1],
                            )
                        )
