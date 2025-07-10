import os
from resistify._loguru import logger


def download_models():
    """
    Manually download models - by default, models will be downloaded automatically when required.
    This command is only necessary if you want to download the models in advance or if you are running in an environment without internet access.
    """
    from transformers import T5EncoderModel, T5Tokenizer
    import esm

    hf_home = os.environ.get("HF_HOME", "~/.cache/huggingface/")
    torch_home = os.environ.get("TORCH_HOME", "~/.cache/torch/")

    logger.info(f"Hugging Face models will be stored under {hf_home}")
    logger.info(f"Torch models will be stored under {torch_home}")

    # ProtT5
    logger.info("Loading ProtT5 models...")
    T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")

    # ESM
    logger.info("Loading ESM models...")
    esm.pretrained.esm2_t33_650M_UR50D()

    logger.info("Models have been downloaded successfully")
