import logging
import esm
from huggingface_hub import snapshot_download

logger = logging.getLogger(__name__)


def download_models():
    logger.info("Downloading ProtT5 (Rostlab/prot_t5_xl_half_uniref50-enc)...")
    snapshot_download("Rostlab/prot_t5_xl_half_uniref50-enc")
    logger.info("ProtT5 downloaded.")

    logger.info("Downloading ESM2-33M (esm2_t33_650M_UR50D via PyTorch Hub)...")
    esm.pretrained.esm2_t33_650M_UR50D()
    logger.info("ESM2-33M downloaded.")

    logger.info("All models downloaded successfully.")
