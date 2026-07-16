import logging
import torch


logger = logging.getLogger(__name__)


def get_device() -> str:
    device = ""
    if torch.cuda.is_available():
        device = "cuda"
    elif torch.backends.mps.is_available():
        device = "mps"
    else:
        device = "cpu"

    logging.info(f"Device selected: {device}")
    return device
