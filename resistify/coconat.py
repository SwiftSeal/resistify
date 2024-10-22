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

log = logging.getLogger(__name__)

class MeanModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = nn.Sequential(
            nn.Linear(4608, 128),
            nn.Dropout(p=0.1),
            nn.ReLU(),
            nn.Linear(128, 4)
        )
        self.softmax = nn.Softmax(dim=1)
    def forward(self, x):
        x = self.softmax(self.linear(x))
        return x

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

def prot_t5_embedding(sequences, lengths, prot_t5_database):
    model = T5EncoderModel.from_pretrained(prot_t5_database)
    tokenizer = T5Tokenizer.from_pretrained(prot_t5_database)
    sequences = [" ".join(re.sub(r"[UZOB]", "X", sequence)) for sequence in sequences]
    ids = tokenizer.batch_encode_plus(sequences, add_special_tokens=True, padding="longest")
    input_ids = torch.tensor(ids['input_ids'])
    attention_mask = torch.tensor(ids['attention_mask'])
    with torch.no_grad():
        embedding_repr = model(input_ids=input_ids, attention_mask=attention_mask)
    
    embeddings = [
        embedding_repr.last_hidden_state[i, :lengths[i]].detach().cpu().numpy()
        for i in range(len(sequences))
    ]
    
    return embeddings

def esm_embedding(sequences, sequence_ids, esm2_database):
    model, alphabet = esm.pretrained.load_model_and_alphabet(esm2_database)
    model.eval()
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
    
    return embeddings

def predict_register_probability(sequences, lengths, merged, registers_model, temp_dir):
    output_path = os.path.join(temp_dir.name, "registers.tsv")
    checkpoint = torch.load(registers_model)
    model = MMModelLSTM()
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    prediction = model(merged, lengths).detach().cpu().numpy()
    with open(output_path, 'w') as outfile:
        for i in range(prediction.shape[0]):
            for j in range(len(sequences[i])):
                prediction_values = " ".join([str(x) for x in prediction[i,j]])
                outfile.write(f"{prediction_values} i\n")
            outfile.write("\n")
    return output_path

def crf(register_path, biocrf_path, crf_model, temp_dir):
    output_path = os.path.join(temp_dir.name, "crf.output.tsv")
    prefix_path = os.path.join(temp_dir.name, "crf.posterior.")
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
            f"{output_path}",
            "-q",
            f"{prefix_path}",
            f"{register_path}",
        ]
    )
    labels, probs = [], []
    current_label = ""
    i = 0
    with open(output_path) as crfo:
        for line in crfo:
            line = line.split()
            if line:
                current_label += line[1]
            else:
                labels.append(current_label)
                probs.append(np.loadtxt(f"{prefix_path.name}_{i}"))
                current_label = ""
                i += 1
    return labels, probs

def coconat(sequences, database, temp_dir, data_dir):
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
    registers_model = os.path.join(data_dir, "dlModel.ckpt")
    biocrf_path = os.path.join(data_dir, bin, "biocrf-static")
    crf_model = os.path.join(data_dir, "crfModel")
    
    log.debug("Extracting N-terminal sequences...")
    sequence_ids, nterminal_sequences, lengths = [], [], []
    for sequence in sequences:
        if sequence[sequence].classification in ["N", "CN"]:
            for annotation in sequence.annotations:
                if annotation.domain == "NB-ARC":
                    nbarc_start = annotation.start
                    break
            nterminal_sequence = sequence[sequence].sequence[:nbarc_start]

            if len(nterminal_sequence) >= 1022:
                log.warning(f"Sequence {sequence[sequence].id} N-terminus too long, skipping. MUST FIX THIS!")
            else:
                sequence_ids.append(sequence)
                nterminal_sequences.append(nterminal_sequence)
                lengths.append(len(nterminal_sequence))
    
    prot_t5_embeddings = prot_t5_embedding(nterminal_sequences, lengths, prot_t5_database)
    esm_embeddings = esm_embedding(nterminal_sequences, sequence_ids, esm2_database)

    merged = [
        torch.from_numpy(np.hstack((prot_t5_embeddings[i], esm_embeddings[i])))
        for i in range(len(nterminal_sequences))
    ]

    merged = torch.nn.utils.rnn.pad_sequence(merged, batch_first=True)

    register_path = predict_register_probability(nterminal_sequences, lengths, merged, registers_model, temp_dir)

    merged = merged.detach().cpu().numpy()

    labels, probs = crf(register_path, biocrf_path, crf_model, temp_dir)

            
    



    