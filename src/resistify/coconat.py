import torch
import esm
from transformers import T5EncoderModel, T5Tokenizer
import regex as re
import numpy as np
import torch.nn as nn
import logging
import os
import warnings
from pathlib import Path
from resistify.annotation import Protein

logger = logging.getLogger(__name__)

warnings.filterwarnings("ignore", category=FutureWarning)


class TransposeX(nn.Module):
    def __init__(self):
        super(TransposeX, self).__init__()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.transpose(x, -2, -1)


# MMModelLSTM and EmbeddingProcessor remain unchanged (except for minor type adjustments)
class MMModelLSTM(nn.Module):
    # ... (MMModelLSTM implementation as before) ...
    def __init__(
        self,
        IN_SIZE: int = 2304,
        OUT_CHANNELS: int = 40,
        KERNEL_SIZE: int = 15,
        LSTM_HIDDEN: int = 128,
        HIDDEN_DIM: int = 64,
        NUM_LSTM: int = 1,
        DROPOUT: float = 0.25,
        OUT_SIZE: int = 8,
    ):
        super().__init__()
        self.cnn = nn.Sequential(
            TransposeX(),
            nn.Conv1d(IN_SIZE, OUT_CHANNELS, KERNEL_SIZE, padding="same"),
            TransposeX(),
            nn.Dropout(p=DROPOUT),
        )
        self.lstm = nn.LSTM(
            OUT_CHANNELS, LSTM_HIDDEN, num_layers=NUM_LSTM, batch_first=True
        )
        self.linear = nn.Sequential(
            nn.Linear(LSTM_HIDDEN, HIDDEN_DIM),
            nn.Dropout(p=DROPOUT),
            nn.ReLU(),
            TransposeX(),
            nn.BatchNorm1d(HIDDEN_DIM),
            TransposeX(),
            nn.Linear(HIDDEN_DIM, OUT_SIZE),
        )
        self.final = nn.Softmax(dim=-1)

    def forward(self, x: torch.Tensor, lengths: list[int]) -> torch.Tensor:
        x = self.cnn(x)
        # Type casting `lengths` to torch.Tensor for pack_padded_sequence
        x = torch.nn.utils.rnn.pack_padded_sequence(
            x, torch.tensor(lengths), batch_first=True, enforce_sorted=False
        )
        x, _ = self.lstm(x)
        x, _ = torch.nn.utils.rnn.pad_packed_sequence(x, batch_first=True)
        x = self.final(self.linear(x))
        return x


class EmbeddingProcessor:
    # ... (EmbeddingProcessor implementation as before) ...
    def __init__(self, device: str):
        self.device = device
        # ProtT5 Model
        self.prot_t5_model: T5EncoderModel = T5EncoderModel.from_pretrained(
            "Rostlab/prot_t5_xl_half_uniref50-enc"
        ).to(self.device)
        self.prot_t5_tokenizer: T5Tokenizer = T5Tokenizer.from_pretrained(
            "Rostlab/prot_t5_xl_half_uniref50-enc"
        )

        # ESM Model
        self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.esm_model.eval()
        self.batch_converter = self.esm_alphabet.get_batch_converter()

    def process_prot_t5_embedding(self, sequence: str, length: int) -> np.ndarray:
        # ... (ProtT5 embedding logic as before) ...
        seq = [" ".join(list(re.sub(r"[UZOB]", "X", sequence)))]

        ids = self.prot_t5_tokenizer.batch_encode_plus(
            seq, add_special_tokens=True, padding="longest"
        )
        input_ids = torch.tensor(ids["input_ids"]).to(self.device)
        attention_mask = torch.tensor(ids["attention_mask"]).to(self.device)

        with torch.no_grad():
            embedding_repr = self.prot_t5_model(
                input_ids=input_ids, attention_mask=attention_mask
            )

        return embedding_repr.last_hidden_state[0, :length].detach().cpu().numpy()

    def process_esm_embedding(self, seq_id: str, seq: str) -> np.ndarray:
        # ... (ESM embedding logic as before) ...
        batch_labels, batch_strs, batch_tokens = self.batch_converter(
            [
                (seq_id, seq),
            ]
        )

        with torch.no_grad():
            results = self.esm_model(
                batch_tokens.to(self.device), repr_layers=[33], return_contacts=False
            )

        token_representations = results["representations"][33]
        seq_len = len(seq)

        return token_representations[0, 1 : seq_len + 1].detach().cpu().numpy()


class CoCoNatPredictor:
    def __init__(self, device: str):
        self.device = torch.device(device)
        logger.info(f"Initializing CoCoNat on device: {self.device}")

        self._load_models()

    def _load_models(self):
        self.embedding_processor = EmbeddingProcessor(device=self.device)

        registers_model = Path(os.path.dirname(__file__)) / "models" / "dlModel.ckpt"
        logger.debug("Loading registers model...")

        checkpoint = torch.load(registers_model, map_location=self.device)

        self.register_model = MMModelLSTM()
        self.register_model.load_state_dict(checkpoint["state_dict"])
        self.register_model.to(self.device)
        self.register_model.eval()
        logger.info("CoCoNat models loaded successfully.")

    def predict_sequence(self, sequence: str, seq_id: str) -> list[float]:
        nterminal_len = len(sequence)

        prot_t5_embeddings = self.embedding_processor.process_prot_t5_embedding(
            sequence, nterminal_len
        )
        esm_embeddings = self.embedding_processor.process_esm_embedding(
            seq_id, sequence
        )

        merged_np = np.hstack((prot_t5_embeddings, esm_embeddings))
        merged_tensor = torch.from_numpy(merged_np).unsqueeze(0).to(self.device).float()

        with torch.no_grad():
            prediction_output = (
                self.register_model(merged_tensor, [nterminal_len])
                .detach()
                .cpu()
                .numpy()
            )

        non_cc_probability = prediction_output[0, :, 0]

        cc_probability = 1 - non_cc_probability

        return cc_probability.tolist()


def predict_coils(proteins: dict[str, Protein], device: str) -> dict[str, Protein]:
    predictor = CoCoNatPredictor(device=device)

    for key, protein in proteins.items():
        logger.debug(f"Processing {protein.id}...")

        nbarc_domains = protein.get_annotation_by_name("NBARC")
        nterminal_start = min([domain.start for domain in nbarc_domains], default=None)

        # Use the N-terminal sequence for prediction
        nterminal_seq = protein.sequence[0 : nterminal_start - 1]

        if nterminal_seq is None:
            logger.debug(f"{protein.id} has no N-terminal sequence, skipping...")
            continue

        nterminal_len = len(nterminal_seq)

        if nterminal_len < 5:
            logger.debug(f"{protein.id} sequence too short for CoCoNat (< 5 residues)")
            continue
        elif nterminal_len >= 1022:
            logger.warning(
                f"{protein.id} sequence quite long (>= 1022), errors might occur."
            )

        try:
            cc_probs_list: list[float] = predictor.predict_sequence(
                sequence=nterminal_seq, seq_id=protein.id
            )

            protein.cc_probs = cc_probs_list

        except Exception as e:
            logger.error(f"Prediction failed for {protein.id}: {e}")

        protein.annotate_cc()

    logger.info("CoCoNat prediction complete.")
    return proteins
