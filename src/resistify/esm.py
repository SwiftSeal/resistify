import logging
import warnings

import torch
from transformers import AutoModel

logger = logging.getLogger(__name__)

logging.getLogger("transformers").setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore", category=FutureWarning)

ESM_MODEL = "Synthyra/ESM2-8M"
ESM_REVISION = "f3c6441"
ESM_DIM = 320


class ESM2Encoder:
    def __init__(self, device: str):
        self.device = torch.device(device)
        self.aa_map = str.maketrans("BJOUZ", "XXXXX")

        logger.info(f"Loading {ESM_MODEL}...")
        model = (
            AutoModel.from_pretrained(
                ESM_MODEL,
                trust_remote_code=True,
                revision=ESM_REVISION,
            )
            .eval()
            .to(self.device)
        )
        self.model = model
        self.tokenizer = model.tokenizer

    @torch.inference_mode()
    def embed(self, sequences: list[str]):
        processed_seqs = [s.upper().translate(self.aa_map) for s in sequences]

        encoded = self.tokenizer(
            processed_seqs, padding="longest", return_tensors="pt"
        ).to(self.device)

        output = self.model(**encoded)

        # Strip BOS token (position 0); remaining positions are residues + EOS [+ PAD]
        embeddings = output.last_hidden_state[:, 1:, :].detach()
        attention_mask = encoded["attention_mask"][:, 1:]

        return embeddings, attention_mask
