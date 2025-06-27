import os
from transformers import T5EncoderModel, T5Tokenizer
import esm
from resistify._loguru import logger


def download_models():
    """
    Downloads the ESM, ProtT5, and ProtTrans models into the specified directory.
    Actually just loads the models, which will download them if not already present.
    """

    hf_home = os.environ.get("HF_HOME", "~/.cache/huggingface/")
    torch_home = os.environ.get("TORCH_HOME", "~/.cache/torch/")

    logger.info(f"Hugging Face models will be stored under {hf_home}")
    logger.info(f"Torch models will be stored under {torch_home}")

    # ProtT5
    logger.info("Loading ProtT5 models...")
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")

    # ESM

    logger.info("Loading ESM models...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()

    logger.info("Models have been downloaded successfully")
