import warnings
import numpy as np
import pickle
import os
import logging
import multiprocessing
import pyhmmer
import threadpoolctl
from concurrent.futures import ProcessPoolExecutor
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

# Number of proteins classified per predict_proba() call/worker task. Batching
# amortizes both the per-call sklearn overhead and process-pool dispatch overhead -
# see _predict_motifs_batch.
BATCH_SIZE = 50

_worker_models: dict[str, object] | None = None


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


def _init_worker(search_type: str):
    # Load the sklearn models once per worker process rather than once per task,
    # and cap BLAS threads to avoid oversubscription across worker processes.
    global _worker_models
    with threadpoolctl.threadpool_limits(limits=1):
        _worker_models = _load_models(search_type)


def _predict_motifs_batch(
    batch: list[tuple[str, np.ndarray]],
) -> list[tuple[str, list[Annotation]]]:
    """
    Classify motifs for a batch of proteins in one process, calling predict_proba
    once per motif for the whole batch (not once per motif per protein) to amortize
    sklearn's per-call overhead across proteins.
    """
    models = _worker_models
    results: dict[str, list[Annotation]] = {seq_id: [] for seq_id, _ in batch}

    for predictor, model in models.items():
        motif_size = MOTIF_SPAN_LENGTHS[predictor]
        window_size = motif_size + 11

        matrices = []
        owners = []
        for seq_id, mat_combined in batch:
            if mat_combined.shape[0] < window_size:
                continue
            windows = np.lib.stride_tricks.sliding_window_view(
                mat_combined, window_size, axis=0
            )
            matrix = windows.transpose(0, 2, 1).reshape(windows.shape[0], -1)
            matrices.append(matrix)
            owners.append((seq_id, matrix.shape[0]))

        if not matrices:
            continue

        proba = model.predict_proba(np.concatenate(matrices, axis=0))

        offset = 0
        for seq_id, n_windows in owners:
            for k, prob in enumerate(proba[offset : offset + n_windows]):
                value = round(float(prob[1]), 4)
                if value > 0.8:
                    i = k + 5
                    results[seq_id].append(
                        Annotation(
                            name=predictor,
                            start=i + 1,
                            end=i + motif_size,
                            type="motif",
                            source="nlrexpress",
                            score=value,
                        )
                    )
            offset += n_windows

    return list(results.items())


def nlrexpress(proteins: dict[str, Protein], threads: int, search_type: str = "all"):
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

    n_batches = -(-len(proteins) // BATCH_SIZE)  # ceil division
    # Split the CPU budget between jackhmmer's own search threads and the motif
    # classification worker processes, mirroring the ~50/50 split that resistify
    # 1.3.0 used between jackhmmer subprocesses and their process pool.
    search_cpus = max(1, threads // 2)
    predict_workers = max(1, min(threads - search_cpus, n_batches))

    batch: list[tuple[str, np.ndarray]] = []
    futures = []

    with ProcessPoolExecutor(
        max_workers=predict_workers,
        initializer=_init_worker,
        initargs=(search_type,),
        mp_context=multiprocessing.get_context("spawn"),
    ) as executor:
        # threadpool limit numpy to 1 to avoid overallocation; motif classification
        # (the CPU-heavy part) now runs in the worker processes, not here.
        with threadpoolctl.threadpool_limits(limits=1):
            for iteration_results in pyhmmer.hmmer.jackhmmer(
                queries,
                database,
                max_iterations=2,
                checkpoints=True,
                cpus=search_cpus,
                E=1e-5,
                domE=1e-5,
            ):
                result1 = iteration_results[0]
                seq_id = result1.hmm.name
                if isinstance(seq_id, bytes):
                    seq_id = seq_id.decode()

                mat1 = -np.log(np.array(result1.hmm.match_emissions) + 1e-8)

                if len(iteration_results) > 1:
                    mat2 = -np.log(
                        np.array(iteration_results[1].hmm.match_emissions) + 1e-8
                    )
                else:
                    mat2 = mat1

                mat_combined = np.concatenate([mat1[1:], mat2[1:]], axis=1)
                batch.append((seq_id, mat_combined))

                if len(batch) >= BATCH_SIZE:
                    futures.append(executor.submit(_predict_motifs_batch, batch))
                    batch = []

            if batch:
                futures.append(executor.submit(_predict_motifs_batch, batch))

        for future in futures:
            for seq_id, annotations in future.result():
                for annotation in annotations:
                    proteins[seq_id].add_annotation(annotation)

    return proteins
