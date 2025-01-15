import torch
import esm
from transformers import T5EncoderModel, T5Tokenizer
import regex as re
import numpy as np
import torch.nn as nn
import subprocess
import logging
import os
import tempfile
import warnings
from resistify._loguru import logger
from resistify.utility import log_percentage

logging.getLogger("transformers").setLevel(logging.CRITICAL)

warnings.filterwarnings("ignore", category=FutureWarning)


class TransposeX(nn.Module):
    def __init__(self):
        super(TransposeX, self).__init__()

    def forward(self, x):
        return torch.transpose(x, -2, -1)


class MMModelLSTM(nn.Module):
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

    def forward(self, x, lengths):
        x = self.cnn(x)
        x = torch.nn.utils.rnn.pack_padded_sequence(
            x, lengths, batch_first=True, enforce_sorted=False
        )
        x, _ = self.lstm(x)
        x, _ = torch.nn.utils.rnn.pad_packed_sequence(x, batch_first=True)
        x = self.final(self.linear(x))
        return x


class EmbeddingProcessor:
    def __init__(self, models_path):
        # Initialize devices and models only once
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if self.device == torch.device("cpu"):
            logger.warning(
                "GPU not available or detected, running on CPU. This will be slower..."
            )

        if models_path is not None:
            # ProtT5 Model
            self.prot_t5_model = T5EncoderModel.from_pretrained(
                os.path.join(models_path, "prott5")
            ).to(self.device)
            self.prot_t5_tokenizer = T5Tokenizer.from_pretrained(
                os.path.join(models_path, "prott5")
            )

            # ESM Model
            self.esm_model, self.esm_alphabet = esm.pretrained.load_model_and_alphabet(
                os.path.join(models_path, "esm", "esm2_t33_650M_UR50D.pt")
            )
        else:
            # ProtT5 Model
            self.prot_t5_model = T5EncoderModel.from_pretrained(
                "Rostlab/prot_t5_xl_half_uniref50-enc"
            ).to(self.device)
            self.prot_t5_tokenizer = T5Tokenizer.from_pretrained(
                "Rostlab/prot_t5_xl_half_uniref50-enc"
            )

            # ESM Model
            self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()

        self.esm_model.eval()
        self.batch_converter = self.esm_alphabet.get_batch_converter()

    def process_prot_t5_embedding(self, sequence, length):
        """
        Compute ProtT5 embeddings for given sequences.

        Args:
            sequences (list): List of protein sequences

        Returns:
            list: List of numpy embeddings
        """
        sequence = [" ".join(list(re.sub(r"[UZOB]", "X", sequence)))]

        ids = self.prot_t5_tokenizer.batch_encode_plus(
            sequence, add_special_tokens=True, padding="longest"
        )
        input_ids = torch.tensor(ids["input_ids"]).to(self.device)
        attention_mask = torch.tensor(ids["attention_mask"]).to(self.device)

        with torch.no_grad():
            embedding_repr = self.prot_t5_model(
                input_ids=input_ids, attention_mask=attention_mask
            )

        embeddings = [
            embedding_repr.last_hidden_state[0, :length].detach().cpu().numpy()
        ]

        return embeddings

    def process_esm_embedding(self, seq_id, seq):
        batch_labels, batch_strs, batch_tokens = self.batch_converter(
            [
                (seq_id, seq),
            ]
        )

        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, repr_layers=[33], return_contacts=False
            )

        token_representations = results["representations"][33]

        embeddings = [
            token_representations[0, 1 : batch_lens[0] - 1].detach().cpu().numpy()
        ]

        return embeddings


def coconat(sequences, models_path: str):
    registers_model = os.path.join(os.path.dirname(__file__), "data", "dlModel.ckpt")
    biocrf_path = os.path.join(os.path.dirname(__file__), "bin", "biocrf-static")
    crf_model = os.path.join(os.path.dirname(__file__), "data", "crfModel")

    logger.debug("Loading embedder...")
    embedding_processor = EmbeddingProcessor(models_path)
    logger.debug("Loading registers model...")
    checkpoint = torch.load(registers_model)
    register_model = MMModelLSTM()
    register_model.load_state_dict(checkpoint["state_dict"])
    register_model.eval()

    temp_dir = tempfile.TemporaryDirectory()
    prediction_file = os.path.join(temp_dir.name, "predictions")
    output_file = os.path.join(temp_dir.name, "out")
    prefix_path = os.path.join(temp_dir.name, "crf")

    total_iterations = len(sequences)
    iteration = 0
    for sequence in sequences:
        logger.debug(f"Processing {sequence.id}...")

        nterminal_seq = sequence.nterminal_sequence

        if nterminal_seq is None:
            logger.debug(f"{sequence.id} has no N-terminus, skipping...")
            continue

        nterminal_len = len(nterminal_seq)

        if nterminal_len < 5:
            logger.debug(f"{sequence.id} N-terminus too short for CoCoNat")
            continue
        elif nterminal_len >= 1022:
            logger.warning(
                f"{sequence.id} N-terminus quite long, errors might occur..."
            )

        prot_t5_embeddings = embedding_processor.process_prot_t5_embedding(
            nterminal_seq, nterminal_len
        )
        esm_embeddings = embedding_processor.process_esm_embedding(
            sequence.id, nterminal_seq
        )

        logger.debug("Merging embeddings")
        merged = [
            torch.from_numpy(np.hstack((prot_t5_embeddings[0], esm_embeddings[0])))
        ]

        merged = torch.nn.utils.rnn.pad_sequence(merged, batch_first=True)

        prediction = register_model(merged, [nterminal_len]).detach().cpu().numpy()

        with open(prediction_file, "w") as outfile:
            for i in range(nterminal_len):
                prediction_values = " ".join([str(x) for x in prediction[0, i]])
                outfile.write(f"{prediction_values} i\n")
            outfile.write("\n")

        logger.debug("Running biocrf")
        subprocess.run(
            [
                f"{biocrf_path}",
                "-test",
                "-m",
                f"{crf_model}",
                "-w",
                "7",
                "-d",
                "posterior-viterbi-sum",
                "-o",
                f"{output_file}",
                "-q",
                f"{prefix_path}",
                f"{prediction_file}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )

        probability_matrix = np.loadtxt(f"{prefix_path}_0")
        cc_probability = 1 - probability_matrix[:, 0]
        sequence.cc_probs = cc_probability
        iteration += 1
        log_percentage(iteration, total_iterations)

    return sequences
