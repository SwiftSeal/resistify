import sys
import re
import shutil
import subprocess
import tempfile
import warnings
import numpy as np
import pickle
import os
import logging
import multiprocessing
import threadpoolctl
import pyhmmer
from concurrent.futures import ProcessPoolExecutor, as_completed
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

# Trailing "-i<N>" suffix that jackhmmer appends to HMM names in checkpoint files
# after the first iteration (e.g. "-i1" after iteration 1).
_ITERATION_SUFFIX = re.compile(r"-i\d+$")


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


_worker_models: dict[str, object] | None = None
_worker_db_path: str | None = None


def _init_worker(search_type: str, db_path: str):
    # Load the sklearn models once per worker process rather than once per task.
    # Cap BLAS threads for the lifetime of this worker process (not just while
    # loading models) - otherwise every predict_proba() call defaults to
    # spinning up one BLAS thread per logical CPU, and with `threads` worker
    # processes doing that concurrently you get threads^2 OS threads fighting
    # over the machine's real cores.
    global _worker_models, _worker_db_path
    threadpoolctl.threadpool_limits(limits=1)
    _worker_models = _load_models(search_type)
    _worker_db_path = db_path


def _read_hmm_checkpoint(path: str) -> dict[str, "pyhmmer.plan7.HMM"]:
    """
    Read a jackhmmer --chkhmm checkpoint file (a plain concatenation of HMMER3
    HMMs, one per query) using pyhmmer's own parser rather than a hand-rolled
    text parser, keyed by query name with any "-i<N>" iteration suffix stripped.
    """
    hmms = {}
    if not os.path.exists(path):
        return hmms
    with pyhmmer.plan7.HMMFile(path) as f:
        for hmm in f:
            name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
            name = _ITERATION_SUFFIX.sub("", name)
            hmms[name] = hmm
    return hmms


def _predict_motifs_batch(
    batch: list[tuple[str, np.ndarray]],
) -> list[tuple[str, list[Annotation]]]:
    """
    Classify motifs for a batch of proteins, calling predict_proba() once per
    motif for the whole batch (not once per motif per protein) to amortize
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


def _search_and_classify_chunk(
    chunk: list[tuple[str, str]],
) -> list[tuple[str, list[Annotation]]]:
    """
    Run jackhmmer as a subprocess against a chunk of query sequences in one
    invocation, read back the resulting checkpoint HMMs, and classify motifs
    for the whole chunk. Runs entirely in a worker process, so it naturally
    overlaps with other chunks being searched/classified in sibling workers.
    """
    seq_ids = [seq_id for seq_id, _ in chunk]

    with tempfile.TemporaryDirectory() as tmp_dir:
        fasta_path = os.path.join(tmp_dir, "query.fa")
        with open(fasta_path, "w") as f:
            for seq_id, sequence in chunk:
                f.write(f">{seq_id}\n{sequence}\n")

        chk_prefix = os.path.join(tmp_dir, "chk")
        cmd = [
            "jackhmmer",
            "--noali",
            "-N",
            "2",
            "--cpu",
            "1",
            "-E",
            "1e-5",
            "--domE",
            "1e-5",
            "--chkhmm",
            chk_prefix,
            fasta_path,
            _worker_db_path,
        ]

        try:
            proc = subprocess.run(cmd, capture_output=True, text=True)
        except FileNotFoundError:
            raise RuntimeError(
                "The 'jackhmmer' executable could not be found. NLRexpress "
                "requires HMMER (http://hmmer.org) to be installed and on PATH."
            )

        if proc.returncode != 0:
            raise RuntimeError(
                f"jackhmmer failed on a chunk of {len(chunk)} sequences "
                f"(starting with {seq_ids[0]!r}): {proc.stderr.strip()}"
            )

        iteration_1 = _read_hmm_checkpoint(f"{chk_prefix}-1.hmm")
        iteration_2 = _read_hmm_checkpoint(f"{chk_prefix}-2.hmm")

    batch: list[tuple[str, np.ndarray]] = []
    for seq_id in seq_ids:
        hmm1 = iteration_1.get(seq_id)
        if hmm1 is None:
            logger.warning(
                f"jackhmmer produced no checkpoint HMM for {seq_id} - skipping."
            )
            continue
        hmm2 = iteration_2.get(seq_id, hmm1)

        mat1 = -np.log(np.array(hmm1.match_emissions) + 1e-8)
        mat2 = -np.log(np.array(hmm2.match_emissions) + 1e-8)
        mat_combined = np.concatenate([mat1[1:], mat2[1:]], axis=1)
        batch.append((seq_id, mat_combined))

    return _predict_motifs_batch(batch)


def _chunk_size_for(n_proteins: int, threads: int) -> int:
    # Aim for a handful of chunks per worker so the process pool load-balances
    # well, while keeping per-chunk jackhmmer invocations large enough to
    # amortize subprocess startup and database-loading overhead.
    return max(1, min(20, n_proteins // max(1, threads * 8)))


def nlrexpress(proteins: dict[str, Protein], threads: int, search_type: str = "all"):
    if shutil.which("jackhmmer") is None:
        logger.critical(
            "The 'jackhmmer' executable was not found on PATH. NLRexpress "
            "requires HMMER (http://hmmer.org) to be installed."
        )
        sys.exit(1)

    db_path = os.path.join(os.path.dirname(__file__), "data", "nlrexpress.fasta")

    items = [(protein.id, protein.sequence) for protein in proteins.values()]
    chunk_size = _chunk_size_for(len(items), threads)
    chunks = [items[i : i + chunk_size] for i in range(0, len(items), chunk_size)]

    logger.info("Running NLRexpress...")

    progress = tqdm(
        total=len(items), desc="NLRexpress", disable=not sys.stdout.isatty()
    )

    with ProcessPoolExecutor(
        max_workers=threads,
        initializer=_init_worker,
        initargs=(search_type, db_path),
        mp_context=multiprocessing.get_context("spawn"),
    ) as executor:
        futures = {
            executor.submit(_search_and_classify_chunk, chunk): len(chunk)
            for chunk in chunks
        }
        for future in as_completed(futures):
            for seq_id, annotations in future.result():
                for annotation in annotations:
                    proteins[seq_id].add_annotation(annotation)
            progress.update(futures[future])

    progress.close()

    return proteins
