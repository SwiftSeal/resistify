import json
import logging
import numpy as np
import torch
from pathlib import Path
from tqdm.auto import tqdm
import xgboost as xgb
from xgboost import XGBClassifier
from resistify.annotation import Protein, Annotation
from resistify.esm import ESM2Encoder

logger = logging.getLogger(__name__)

MODELS_DIR = Path(__file__).parent / "data" / "models"

WINDOW_SIZE = 30

MOTIF_SPAN_LENGTHS = {
    "VG": 5,
    "P-loop": 9,
    "RNBS-A": 10,
    "RNBS-B": 7,
    "RNBS-C": 10,
    "RNBS-D": 9,
    "Walker-B": 8,
    "GLPL": 5,
    "MHD": 3,
    "extEDVID": 12,
    #    "aA": 7,
    #    "aC": 6,
    #    "aD3": 13,
    #    "bA": 10,
    #    "bC": 8,
    #    "bDaD1": 16,
    "LxxLxL": 6,
}


def _load_models(models_dir: Path, search_type: str, threads: int):
    models = {}
    motifs = list(MOTIF_SPAN_LENGTHS.keys()) if search_type == "all" else [search_type]

    for motif in motifs:
        model_path = models_dir / f"{motif}.ubj"
        meta_path = models_dir / f"{motif}_meta.json"

        if not model_path.exists():
            logger.warning(f"No model found for motif {motif}, skipping")
            continue

        clf = XGBClassifier()
        clf.load_model(model_path)
        clf.set_params(nthread=threads)

        with open(meta_path) as f:
            meta = json.load(f)

        models[motif] = {"model": clf, "meta": meta}
        logger.debug(f"Loaded {motif} (threshold={meta['threshold']:.4f})")

    return models



def _make_windows(matrix: np.ndarray, window_size: int) -> xgb.DMatrix:
    pad = window_size // 2
    padded = np.pad(matrix, ((pad, pad), (0, 0)), mode="constant")
    windowed = np.lib.stride_tricks.sliding_window_view(padded, window_size, axis=0)[
        : len(matrix)
    ]
    windowed = windowed.transpose(0, 2, 1).reshape(len(matrix), -1)
    # convert to dmatrix once for speeeeeed
    return xgb.DMatrix(windowed)


def nlrexpress(
    proteins: dict[str, Protein],
    search_type: str = "all",
    device: str = "cpu",
    threads: int = 1,
    encoder: ESM2Encoder | None = None,
):
    logger.info(f"Running motif classifier for '{search_type}' motifs")

    models = _load_models(MODELS_DIR, search_type, threads)
    if not models:
        logger.warning("No models loaded, skipping")
        return proteins

    torch.set_num_threads(threads)
    if encoder is None:
        encoder = ESM2Encoder(device)

    for seq_id, protein in tqdm(proteins.items(), desc="Predicting motifs"):
        if protein.length < WINDOW_SIZE:
            continue

        embeddings, _ = encoder.embed([protein.sequence])
        emb = embeddings[0, :protein.length, :].cpu().float().numpy()
        if np.isnan(emb).any():
            raise RuntimeError(f"NaN embeddings found for sequence {seq_id}")

        windows = _make_windows(emb, WINDOW_SIZE)

        for motif, clf in models.items():
            threshold = clf["meta"]["threshold"]
            span = MOTIF_SPAN_LENGTHS[motif]

            proba = clf["model"].get_booster().predict(windows)

            for idx in np.where(proba >= threshold)[0]:
                end = int(idx + span)
                if end > protein.length:
                    continue
                protein.add_annotation(
                    Annotation(
                        name=motif,
                        type="motif",
                        start=int(idx + 1),
                        end=end,
                        source="motif_classifier",
                        score=float(proba[idx]),
                    )
                )

    logger.info("Motif classification completed")
    return proteins
