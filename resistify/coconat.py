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

log = logging.getLogger(__name__)
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

    def forward(self, x, l):
        x = self.cnn(x)
        x = torch.nn.utils.rnn.pack_padded_sequence(
            x, l, batch_first=True, enforce_sorted=False
        )
        x, _ = self.lstm(x)
        x, _ = torch.nn.utils.rnn.pad_packed_sequence(x, batch_first=True)
        x = self.final(self.linear(x))
        return x


class EmbeddingProcessor:
    def __init__(self, models_path):
        # Initialize devices and models only once
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        log.debug(f"Device is {self.device}")

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

    def process_prot_t5_embedding(self, sequences):
        """
        Compute ProtT5 embeddings for given sequences.

        Args:
            sequences (list): List of protein sequences

        Returns:
            list: List of numpy embeddings
        """
        lengths = [len(sequence) for sequence in sequences]
        sequences = [
            " ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequences
        ]

        ids = self.prot_t5_tokenizer.batch_encode_plus(
            sequences, add_special_tokens=True, padding="longest"
        )
        input_ids = torch.tensor(ids["input_ids"]).to(self.device)
        attention_mask = torch.tensor(ids["attention_mask"]).to(self.device)

        with torch.no_grad():
            embedding_repr = self.prot_t5_model(
                input_ids=input_ids, attention_mask=attention_mask
            )

        embeddings = [
            embedding_repr.last_hidden_state[i, : lengths[i]].detach().cpu().numpy()
            for i in range(len(sequences))
        ]

        return embeddings

    def process_esm_embedding(self, chunk_ids, chunk_seqs):
        batch_labels, batch_strs, batch_tokens = self.batch_converter(
            list(zip(chunk_ids, chunk_seqs))
        )

        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = self.esm_model(
                batch_tokens, repr_layers=[33], return_contacts=False
            )

        token_representations = results["representations"][33]

        embeddings = [
            token_representations[i, 1 : tokens_len - 1].detach().cpu().numpy()
            for i, tokens_len in enumerate(batch_lens)
        ]

        return embeddings


def coconat(sequences, models_path: str):
    registers_model = os.path.join(os.path.dirname(__file__), "data", "dlModel.ckpt")
    biocrf_path = os.path.join(os.path.dirname(__file__), "bin", "biocrf-static")
    crf_model = os.path.join(os.path.dirname(__file__), "data", "crfModel")

    log.debug("Loading embedder...")
    embedding_processor = EmbeddingProcessor(models_path)
    log.debug("Loading registers model...")
    checkpoint = torch.load(registers_model)
    register_model = MMModelLSTM()
    register_model.load_state_dict(checkpoint["state_dict"])
    register_model.eval()

    for chunk_start in range(0, len(sequences), 5):
        log.debug("Processing batch...")

        chunk_ids, chunk_seqs, chunk_lengths = [], [], []
        for sequence in sequences[chunk_start : chunk_start + 5]:
            n_terminal_seq = sequence.get_nterminal()
            if n_terminal_seq is None:
                log.debug(f"{sequence.id} has no N-terminus, skipping...")
                continue
            elif len(n_terminal_seq) < 5:
                log.debug(f"{sequence.id} N-terminus too short for CoCoNat")
                continue
            elif len(n_terminal_seq) >= 1022:
                log.warning(
                    f"{sequence.id} N-terminus quite long, errors might occur..."
                )

            chunk_ids.append(sequence.id)
            chunk_seqs.append(n_terminal_seq)
            chunk_lengths.append(len(n_terminal_seq))

        if len(chunk_ids) == 0:
            continue

        prot_t5_embeddings = embedding_processor.process_prot_t5_embedding(chunk_seqs)
        esm_embeddings = embedding_processor.process_esm_embedding(
            chunk_ids, chunk_seqs
        )

        merged = [
            torch.from_numpy(np.hstack((prot_t5_embeddings[i], esm_embeddings[i])))
            for i in range(len(chunk_ids))
        ]

        merged = torch.nn.utils.rnn.pad_sequence(merged, batch_first=True)

        prediction = register_model(merged, chunk_lengths).detach().cpu().numpy()

        temp_dir = tempfile.TemporaryDirectory()
        prediction_file = os.path.join(temp_dir.name, "predictions")
        output_file = os.path.join(temp_dir.name, "out")
        prefix_path = os.path.join(temp_dir.name, "crf")

        with open(prediction_file, "w") as outfile:
            for i in range(prediction.shape[0]):
                for j in range(chunk_lengths[i]):
                    prediction_values = " ".join([str(x) for x in prediction[i, j]])
                    outfile.write(f"{prediction_values} i\n")
                outfile.write("\n")

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

        cc_probabilities = {}

        for i, sequence_id in enumerate(chunk_ids):
            log.debug(
                f"Loading crf probabilities for {sequence_id} in {prefix_path}_{i}"
            )
            probability_matrix = np.loadtxt(f"{prefix_path}_{i}")
            # extract first column
            cc_probability = 1 - probability_matrix[:, 0]
            cc_probabilities[sequence_id] = cc_probability

        for sequence in sequences:
            if sequence.id in cc_probabilities.keys():
                sequence.cc_probs = cc_probabilities[sequence.id]

    return sequences
