import torch
import esm
from transformers import T5EncoderModel, T5Tokenizer
import regex as re
import numpy as np
import torch.nn as nn
import subprocess
import logging
import sys
import os
import tempfile
import csv
from resistify.annotations import Annotation
import warnings

log = logging.getLogger(__name__)
logging.getLogger("transformers").setLevel(logging.CRITICAL)


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


def split_sequences(sequences):
    sequence_ids, lengths, chunk_ids, chunk_sequences = [], [], [], []
    for sequence in sequences:
        if sequence.classification in ["N", "CN", "NL", "CNL"]:
            nterminal_sequence = sequence.get_nterminal()
            # Probably arbitrary, but 5 is unlikely to affect results
            if len(nterminal_sequence) < 5:
                log.debug(f"{sequence.id} N-terminal too short, will not be processed")
                continue
            sequence_ids.append(sequence.id)
            lengths.append(len(nterminal_sequence))
            if len(nterminal_sequence) >= 1022:
                log.debug(
                    f"N-terminus of {sequence.id} is too long, splitting into chunks..."
                )
                for i in range(0, len(nterminal_sequence), 1022):
                    chunk_ids.append(f"{sequence.id}_{i}")
                    chunk_sequences.append(nterminal_sequence[i : i + 1022])
            else:
                chunk_ids.append(f"{sequence.id}_0")
                chunk_sequences.append(nterminal_sequence)
    return sequence_ids, lengths, chunk_ids, chunk_sequences


class EmbeddingProcessor:
    def __init__(self):
        # Initialize devices and models only once
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        
        # ProtT5 Model
        self.prot_t5_model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc").to(self.device)
        self.prot_t5_tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
        
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
            embedding_repr = self.prot_t5_model(input_ids=input_ids, attention_mask=attention_mask)

        embeddings = [
            embedding_repr.last_hidden_state[i, : lengths[i]].detach().cpu().numpy()
            for i in range(len(sequences))
        ]
        
        return embeddings

    def process_esm_embedding(self, sequences, sequence_ids):
        """
        Compute ESM embeddings for given sequences.
        
        Args:
            sequences (list): List of protein sequences
            sequence_ids (list): List of sequence identifiers
        
        Returns:
            list: List of numpy embeddings
        """
        batch_labels, batch_strs, batch_tokens = self.batch_converter(
            list(zip(sequence_ids, sequences))
        )

        batch_lens = (batch_tokens != self.esm_alphabet.padding_idx).sum(1)

        with torch.no_grad():
            results = self.esm_model(batch_tokens, repr_layers=[33], return_contacts=False)

        token_representations = results["representations"][33]

        embeddings = [
            token_representations[i, 1 : tokens_len - 1].detach().cpu().numpy()
            for i, tokens_len in enumerate(batch_lens)
        ]
        
        return embeddings

def merge_embeddings(sequence_ids, chunk_ids, prot_t5_embeddings, esm_embeddings):
    """
    Merge embeddings for chunked sequences.
    
    Args:
        sequence_ids (list): List of original sequence IDs
        chunk_ids (list): List of chunk IDs
        prot_t5_embeddings (list): List of ProtT5 embeddings
        esm_embeddings (list): List of ESM embeddings
    
    Returns:
        list: List of merged torch embeddings
    """
    previous = None
    merged_prot, merged_esm = [], []
    for i, chunk_id in enumerate(chunk_ids):
        sequence_id = chunk_id.rsplit("_", 1)[0]
        if sequence_id != previous:
            merged_prot.append(prot_t5_embeddings[i])
            merged_esm.append(esm_embeddings[i])
        else:
            merged_prot[-1] = np.vstack((merged_prot[-1], prot_t5_embeddings[i]))
            merged_esm[-1] = np.vstack((merged_esm[-1], esm_embeddings[i]))
        previous = sequence_id
    
    merged = [
        torch.from_numpy(np.hstack((merged_prot[i], merged_esm[i])))
        for i in range(len(sequence_ids))
    ]
    return merged

def coconat_embeddings(sequences, embedding_processor, register_model):
    """
    Compute embeddings for sequences using pre-loaded models.
    
    Args:
        sequences (list): List of sequences to process
        embedding_processor (EmbeddingProcessor): Initialized embedding processor
    
    Returns:
        tuple: Processed sequences, merged embeddings, lengths, chunk_ids
    """
    biocrf_path = os.path.join(os.path.dirname(__file__), "bin", "biocrf-static")
    crf_model = os.path.join(os.path.dirname(__file__), "data", "crfModel")

    sequence_ids, lengths, chunk_ids, chunk_sequences = split_sequences(sequences)

    # Compute embeddings
    prot_t5_embeddings = embedding_processor.process_prot_t5_embedding(chunk_sequences)
    esm_embeddings = embedding_processor.process_esm_embedding(chunk_sequences, chunk_ids)

    # Merge embeddings
    merged = merge_embeddings(sequence_ids, chunk_ids, prot_t5_embeddings, esm_embeddings)
    merged = torch.nn.utils.rnn.pad_sequence(merged, batch_first=True)

    prediction = register_model(merged, lengths).detach().cpu().numpy()

    temp_dir = tempfile.TemporaryDirectory()
    prediction_file = os.path.join(temp_dir.name, "predictions")
    output_file = os.path.join(temp_dir.name, "out")
    prefix_path = os.path.join(temp_dir.name, "crf")

    with open(prediction_file, "w") as outfile:
        for i in range(prediction.shape[0]):
            for j in range(lengths[i]):
                prediction_values = " ".join([str(x) for x in prediction[i, j]])
                outfile.write(f"{prediction_values} i\n")
            outfile.write("\n")
        
    subprocess.call(
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
        #stdout=subprocess.PIPE,
        #stderr=subprocess.PIPE,
    )

    cc_probabilities = {}

    for i, sequence_id in enumerate(sequence_ids):
        log.debug(f"Loading crf probabilities for {sequence_id} in {prefix_path}_{i}")
        probability_matrix = np.loadtxt(f"{prefix_path}_{i}")
        # extract first column
        cc_probability = 1 - probability_matrix[:, 0]
        cc_probabilities[sequence_id] = cc_probability


    return cc_probabilities

def coconat(sequences):
    registers_model = os.path.join(os.path.dirname(__file__), "data", "dlModel.ckpt")

    log.debug("Loading embedder...")
    embedding_processor = EmbeddingProcessor()
    log.debug("Loading registers model...")
    checkpoint = torch.load(registers_model)
    register_model = MMModelLSTM()
    register_model.load_state_dict(checkpoint["state_dict"])
    register_model.eval()
    
    for i in range(0, len(sequences), 5):
        log.debug("Processing batch...")
        batch = sequences[i:i+5]
        cc_probabilities = coconat_embeddings(batch, embedding_processor, register_model)
        for sequence in sequences:
            if sequence.id in cc_probabilities.keys():
                sequence.cc_probs = cc_probabilities[sequence.id]

    return sequences




