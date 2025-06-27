import os
from transformers import T5EncoderModel, T5Tokenizer
import esm
from resistify._loguru import logger


def download_models(models_directory: str):
    """
    Downloads the ESM, ProtT5, and ProtTrans models into the specified directory.
    Actually just loads the models, which will download them if not already present.
    """
    # Create the directory if it doesn't exist
    logger.info(f"Creating models directory at {models_directory}")
    os.makedirs(models_directory, exist_ok=True)

    # set environment variables to point to the models directory
    os.environ["HF_HOME"] = models_directory
    os.environ["TORCH_HOME"] = models_directory
    
    # ProtT5
    logger.info("Loading ProtT5 models...")
    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")

    # ESM

    logger.info("Loading ESM models...")
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()

    logger.info(f"Models downloaded successfully. \nSupply these with the argument `--models {models_directory}`\n Alternatively, you can set the environment variables HF_HOME and TORCH_HOME to {models_directory}")