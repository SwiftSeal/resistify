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
            OUT_SIZE: int = 8
        ):
        super().__init__()
        self.cnn = nn.Sequential(
            TransposeX(),
            nn.Conv1d(
                IN_SIZE,
                OUT_CHANNELS,
                KERNEL_SIZE,
                padding='same'
            ),
            TransposeX(),
            nn.Dropout(p=DROPOUT)
        )
        self.lstm = nn.LSTM(
            OUT_CHANNELS,
            LSTM_HIDDEN,
            num_layers=NUM_LSTM,
            batch_first=True
        )
        self.linear = nn.Sequential(
            nn.Linear(LSTM_HIDDEN, HIDDEN_DIM),
            nn.Dropout(p=DROPOUT),
            nn.ReLU(),
            TransposeX(),
            nn.BatchNorm1d(HIDDEN_DIM),
            TransposeX(),
            nn.Linear(HIDDEN_DIM, OUT_SIZE)
        )
        self.final = nn.Softmax(dim=-1)
    
    def forward(self, x, l):
        x = self.cnn(x)
        x = torch.nn.utils.rnn.pack_padded_sequence(x, l, batch_first=True, enforce_sorted=False)
        x, _ = self.lstm(x)
        x, _ = torch.nn.utils.rnn.pad_packed_sequence(x, batch_first=True)
        x = self.final(self.linear(x))
        return x

PROT_T5_FILES = [
    "config.json",
    "pytorch_model.bin",
    "special_tokens_map.json",
    "spiece.model",
    "tokenizer_config.json",
]

ESM2_FILES = [
    "esm2_t33_650M_UR50D-contact-regression.pt",
    "esm2_t33_650M_UR50D.pt"
]

def split_sequences(sequences):
    sequence_ids, lengths, chunk_ids, chunk_sequences = [], [], [], []
    for sequence in sequences:
        if sequence.classification in ["N", "CN"]:
            nterminal_sequence = sequence.get_nterminal()
            if len(nterminal_sequence) == 0:
                log.debug(f"{sequence.id} has no N-terminal region!")
                continue
            sequence_ids.append(sequence.id)
            lengths.append(len(nterminal_sequence))
            if len(nterminal_sequence) >= 1022:
                log.debug(f"N-terminus of {sequence.id} is too long, splitting into chunks...")
                for i in range(0, len(nterminal_sequence), 1022):
                    chunk_ids.append(f"{sequence.id}_{i}")
                    chunk_sequences.append(nterminal_sequence[i:i+1022])
            else:
                chunk_ids.append(f"{sequence.id}_0")
                chunk_sequences.append(nterminal_sequence)
    return sequence_ids, lengths, chunk_ids, chunk_sequences

def prot_t5_embedding(sequences, prot_t5_database):
    lengths = [len(sequence) for sequence in sequences]
    log.debug("Loading ProtT5 model...")
    # suppress warning about weights
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        model = T5EncoderModel.from_pretrained(prot_t5_database)
        tokenizer = T5Tokenizer.from_pretrained(prot_t5_database)
    log.debug("Generating ProtT5 embeddings...")
    sequences = [" ".join(list(re.sub(r"[UZOB]", "X", sequence))) for sequence in sequences]
    ids = tokenizer.batch_encode_plus(sequences, add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids'])
    attention_mask = torch.tensor(ids['attention_mask'])
    log.debug("Computing ProtT5 embeddings...")
    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
    
    embeddings = [
        embedding_repr.last_hidden_state[i, :lengths[i]].detach().cpu().numpy()
        for i in range(len(sequences))
    ]
    log.debug(f"Embedding for Prot5 completed, length of embeddings is {len(embeddings)}")
    return embeddings

def esm_embedding(sequences, sequence_ids, esm2_database):
    log.debug("Loading ESM2 model...")
    # suppress warning about weights
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FutureWarning)
        model, alphabet = esm.pretrained.load_model_and_alphabet(esm2_database)
    model.eval()
    log.debug("Generating ESM2 embeddings...")
    batch_converter = alphabet.get_batch_converter()
    batch_labels, batch_strs, batch_tokens = batch_converter(list(zip(sequence_ids, sequences)))
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=False)
    token_representations = results["representations"][33]
    embeddings = [
        token_representations[i, 1 : tokens_len - 1].detach().cpu().numpy()
        for i, tokens_len in enumerate(batch_lens)
    ]
    log.debug(f"Embedding for ESM2 completed, length of embeddings is {len(embeddings)}")
    return embeddings

def merge_embeddings(chunk_ids, embeddings):
    previous = None
    merged_embedding = []
    for i, chunk_id in enumerate(chunk_ids):
        sequence_id = chunk_id.rsplit("_", 1)[0]
        if sequence_id != previous:
            log.debug(f"Creating new merged embedding for {sequence_id}")
            merged_embedding.append(embeddings[i])
        else:
            log.debug(f"Merging embedding for {sequence_id}")
            merged_embedding[-1] = np.vstack((merged_embedding[-1], embeddings[i]))
        previous = sequence_id
    log.debug(f"Merged embedding length: {len(merged_embedding)}")
    return merged_embedding


def predict_register_probability(lengths, merged, registers_model):
    output_path = tempfile.NamedTemporaryFile()
    log.debug("Loading registers model...")
    checkpoint = torch.load(registers_model, weights_only=False)
    model = MMModelLSTM()
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    prediction = model(merged, lengths).detach().cpu().numpy()
    with open(output_path.name, 'w') as outfile:
        for i in range(prediction.shape[0]):
            for j in range(lengths[i]):
                prediction_values = " ".join([str(x) for x in prediction[i,j]])
                outfile.write(f"{prediction_values} i\n")
            outfile.write("\n")
    return output_path

def crf(sequence_ids, register_path, biocrf_path, crf_model):
    output_path = tempfile.NamedTemporaryFile()
    prefix_directory = tempfile.TemporaryDirectory()
    prefix_path = os.path.join(prefix_directory.name, "crf")
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
            f"{output_path.name}",
            "-q",
            f"{prefix_path}",
            f"{register_path.name}",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    # A prefix file is created for each sequence (or sequence chunk).
    # Iterate through the sequence ids and extract the first column of the prefix file - these are the per-residue probabilities.
    cc_probabilities = {}
    
    for i, sequence_id in enumerate(sequence_ids):
        log.debug(f"Loading crf probabilities for {sequence_id} in {prefix_path}_{i}")
        probability_matrix = np.loadtxt(f"{prefix_path}_{i}")
        # extract first column
        cc_probability = probability_matrix[:, 0]
        cc_probabilities[sequence_id] = cc_probability

    return cc_probabilities

def coconat(sequences, database):
    """
    Use Coconat to predict coiled-coil domains in the N-terminal regions of NLRs.
    """
    log.debug("Validating Coconat database...")
    if not os.path.isdir(database):
        log.error(f"Coconat database not found at {database}")
        sys.exit(1)
    
    for file in PROT_T5_FILES:
        if not os.path.isfile(os.path.join(database, "prot_t5_xl_uniref50", file)):
            log.error(f"Coconat database is missing file: {file}")
            sys.exit(1)
    
    for file in ESM2_FILES:
        if not os.path.isfile(os.path.join(database, "esm2", file)):
            log.error(f"Coconat database is missing file: {file}")
            sys.exit(1)
    
    prot_t5_database = os.path.join(database, "prot_t5_xl_uniref50")
    esm2_database = os.path.join(database, "esm2", "esm2_t33_650M_UR50D.pt")
    registers_model = os.path.join(os.path.dirname(__file__), "data", "dlModel.ckpt")
    biocrf_path = os.path.join(os.path.dirname(__file__), "bin", "biocrf-static")
    crf_model = os.path.join(os.path.dirname(__file__), "data", "crfModel")
    
    log.debug("Extracting N-terminal sequences...")
    sequence_ids, lengths, chunk_ids, chunk_sequences = split_sequences(sequences)

    
    prot_t5_embeddings = prot_t5_embedding(chunk_sequences, prot_t5_database)
    prot_t5_embeddings = merge_embeddings(chunk_ids, prot_t5_embeddings)

    esm_embeddings = esm_embedding(chunk_sequences, chunk_ids, esm2_database)
    esm_embeddings = merge_embeddings(chunk_ids, esm_embeddings)

    merged = [
        torch.from_numpy(np.hstack((prot_t5_embeddings[i], esm_embeddings[i])))
        for i in range(len(sequence_ids))
    ]

    merged = torch.nn.utils.rnn.pad_sequence(merged, batch_first=True)

    register_path = predict_register_probability(lengths, merged, registers_model)

    merged = merged.detach().cpu().numpy()

    cc_probabilities = crf(sequence_ids, register_path, biocrf_path, crf_model)

    #with open("coconat_results.txt", "w") as outfile:
    #    result_writer = csv.writer(outfile, delimiter="\t")
    #    result_writer.writerow(["sequence_id", "position", "cc_probability"])
    #    for i in range(len(nterminal_sequences)):
    #        for j in range(len(nterminal_sequences[i])):
    #            result_writer.writerow(
    #                [
    #                    sequence_ids[i],
    #                    j,
    #                    probs[i][j][0],
    #                ]
    #            )
    
    for sequence in sequences:
        if sequence.id in cc_probabilities:
            sequence.cc_probs = cc_probabilities[sequence.id]

    return sequences

            
    



    
