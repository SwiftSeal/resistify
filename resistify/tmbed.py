# Copyright 2022 Rostlab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import os
import numpy
import math
import random
import torch
import torch.nn as nn
import torch.nn.functional as F

import logging

from transformers import T5EncoderModel, T5Tokenizer

log = logging.getLogger(__name__)
logging.getLogger("transformers").setLevel(logging.CRITICAL)


class Decoder(nn.Module):
    def __init__(self):
        super().__init__()

        self._init_transitions()

    def _init_transitions(self):
        num_tags = 27

        end_transitions = torch.full((num_tags,), -100)
        start_transitions = torch.full((num_tags,), -100)

        transitions = torch.full((num_tags, num_tags), -100)

        for i in [0, 5, 10, 15, 20, -2, -1]:
            start_transitions[i] = 0  # B1a, B1b, H1a, H1b, S1, i, o

        for i in range(4):
            transitions[0 + i, 1 + i] = 0  # Bxa -> Bya
            transitions[5 + i, 6 + i] = 0  # Bxb -> Byb
            transitions[10 + i, 11 + i] = 0  # Hxa -> Hya
            transitions[15 + i, 16 + i] = 0  # Hxb -> Hyb
            transitions[20 + i, 21 + i] = 0  # Sx  -> Sy

        for i in [4, 9, 14, 19, 24]:
            transitions[i, i] = 0  # X5 -> X5

        transitions[4, -1] = 0  # B5a -> o
        transitions[9, -2] = 0  # B5b -> i
        transitions[14, -1] = 0  # H5a -> o
        transitions[19, -2] = 0  # H5b -> i
        transitions[24, -2:] = 0  # S5  -> (i, o)

        transitions[-2, 0] = 0  # i -> B1a
        transitions[-2, 10] = 0  # i -> H1a
        transitions[-2, -2:] = 0  # i -> (i, o)

        transitions[-1, 5] = 0  # o -> B1b
        transitions[-1, 15] = 0  # o -> H1b
        transitions[-1, -2:] = 0  # o -> (i, o)

        for i in [4, 9, 14, 19, 24, -2, -1]:
            end_transitions[i] = 0  # B5a, B5b, H5a, H5b, S5, i, o

        repeats = torch.tensor([10, 10, 5, 1, 1], dtype=torch.int32)

        mapping = torch.arange(7, dtype=torch.int32)
        mapping = mapping.repeat_interleave(
            torch.tensor([5, 5, 5, 5, 5, 1, 1])  # B  # H  # S  # i
        )  # o

        assert repeats.sum() == num_tags
        assert mapping.shape == (num_tags,)

        self.register_buffer("transitions", tensor=transitions)
        self.register_buffer("end_transitions", tensor=end_transitions)
        self.register_buffer("start_transitions", tensor=start_transitions)

        self.register_buffer("repeats", tensor=repeats)
        self.register_buffer("mapping", tensor=mapping)

    def forward(self, emissions, mask):
        mask = mask.transpose(0, 1).bool()

        emissions = emissions.permute(2, 0, 1)
        emissions = emissions.repeat_interleave(self.repeats, dim=2)

        decoded = self._viterbi_decode(emissions, mask)
        decoded = self.mapping[decoded]

        return decoded

    def _viterbi_decode(self, emissions, mask):
        device = emissions.device

        seq_length, batch_size, num_tags = emissions.shape

        score = self.start_transitions + emissions[0]

        history = torch.zeros(
            (seq_length, batch_size, num_tags), dtype=torch.long, device=device
        )

        for i in range(1, seq_length):
            next_score = (
                self.transitions + score.unsqueeze(2) + emissions[i].unsqueeze(1)
            )

            next_score, indices = next_score.max(dim=1)

            score = torch.where(mask[i].unsqueeze(-1), next_score, score)

            history[i - 1] = indices

        score = score + self.end_transitions

        _, end_tag = score.max(dim=1)

        seq_ends = mask.long().sum(dim=0) - 1

        history = history.transpose(1, 0)

        history.scatter_(
            1,
            seq_ends.view(-1, 1, 1).expand(-1, 1, num_tags),
            end_tag.view(-1, 1, 1).expand(-1, 1, num_tags),
        )

        history = history.transpose(1, 0)

        best_tags = torch.zeros((batch_size, 1), dtype=torch.long, device=device)

        best_tags_arr = torch.zeros(
            (seq_length, batch_size), dtype=torch.long, device=device
        )

        for idx in range(seq_length - 1, -1, -1):
            best_tags = torch.gather(history[idx], 1, best_tags)

            best_tags_arr[idx] = best_tags.view(batch_size)

        return best_tags_arr.transpose(0, 1)


class T5Encoder:
    def __init__(self, models_path, use_gpu=True):
        if use_gpu and torch.cuda.is_available():
            self._load_models(models_path, torch.float16)
            self.encoder_model = self.encoder_model.eval().cuda()
        else:
            self._load_models(models_path, torch.float32)
            self.encoder_model = self.encoder_model.eval()

        self.aa_map = str.maketrans("BJOUZ", "XXXXX")

    def _load_models(self, models_path, dtype):
        if models_path is not None:
            self.tokenizer = T5Tokenizer.from_pretrained(
                os.path.join(models_path, "prott5"),
                do_lower_case=False,
            )
            self.encoder_model = T5EncoderModel.from_pretrained(
                os.path.join(models_path, "prott5"),
                torch_dtype=dtype,
            )
        else:
            self.tokenizer = T5Tokenizer.from_pretrained(
                "Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False
            )
            self.encoder_model = T5EncoderModel.from_pretrained(
                "Rostlab/prot_t5_xl_half_uniref50-enc", torch_dtype=dtype
            )

    def device(self):
        return self.encoder_model.device

    def to_cpu(self):
        self.encoder_model = self.encoder_model.cpu().float()

    def to_cuda(self):
        self.encoder_model = self.encoder_model.half().cuda()

    def embed(self, sequence):
        sequence = sequence.upper().translate(self.aa_map)

        tokens = [" ".join(list(sequence))]
        tokens = self.tokenizer.batch_encode_plus(
            tokens, padding="longest", add_special_tokens=True
        )

        device = self.encoder_model.device
        input_ids = torch.tensor(tokens["input_ids"], device=device)
        attention_mask = torch.tensor(tokens["attention_mask"], device=device)

        with torch.no_grad():
            embeddings = self.encoder_model(
                input_ids=input_ids, attention_mask=attention_mask
            )

        return embeddings.last_hidden_state.detach()


class SeqNorm(nn.Module):
    def __init__(self, channels, eps=1e-6, affine=True):
        super().__init__()

        if affine:
            self.bias = nn.Parameter(torch.zeros(1, channels, 1, 1))
            self.weight = nn.Parameter(torch.ones(1, channels, 1, 1))
        else:
            self.register_parameter("bias", None)
            self.register_parameter("weight", None)

        self.register_buffer(name="eps", tensor=torch.tensor(float(eps)))
        self.register_buffer(name="channels", tensor=torch.tensor(channels))

    def forward(self, x, mask):
        mask_rsum = 1.0 / (mask.sum(dim=(2, 3), keepdims=True) * self.channels)

        x = x * mask

        mean = x.sum(dim=(1, 2, 3), keepdims=True) * mask_rsum

        x = (x - mean) * mask

        var = x.square().sum(dim=(1, 2, 3), keepdims=True) * mask_rsum

        x = x * torch.rsqrt(var + self.eps)

        if self.weight is not None:
            x = (x * self.weight) + self.bias

        return x


class Conv(nn.Module):
    def __init__(
        self,
        in_channels,
        out_channels,
        kernel_size,
        stride=1,
        padding=0,
        dilation=1,
        groups=1,
        padding_mode="zeros",
    ):
        super().__init__()

        self.func = nn.ReLU(inplace=True)

        self.norm = SeqNorm(channels=out_channels, eps=1e-6, affine=True)

        self.conv = nn.Conv2d(
            in_channels=in_channels,
            out_channels=out_channels,
            kernel_size=(kernel_size, 1),
            stride=(stride, 1),
            padding=(padding, 0),
            dilation=(dilation, 1),
            groups=groups,
            bias=False,
            padding_mode=padding_mode,
        )

        # Init Conv Params
        nn.init.xavier_uniform_(self.conv.weight)

    def forward(self, x, mask):
        x = self.func(self.norm(self.conv(x), mask))

        return x * mask


class CNN(nn.Module):
    def __init__(self, channels):
        super().__init__()

        self.input = Conv(1024, channels, 1, 1, 0)

        self.dwc1 = Conv(channels, channels, 9, 1, 4, groups=channels)
        self.dwc2 = Conv(channels, channels, 21, 1, 10, groups=channels)

        self.dropout = nn.Dropout2d(p=0.50, inplace=True)

        self.output = nn.Conv2d(3 * channels, 5, 1, 1, 0)

        # Init Output Params
        nn.init.zeros_(self.output.bias)
        nn.init.xavier_uniform_(self.output.weight)

    def forward(self, x, mask):
        x = self.input(x, mask)

        z1 = self.dwc1(x, mask)
        z2 = self.dwc2(x, mask)

        x = torch.cat([x, z1, z2], dim=1)

        x = self.dropout(x)

        x = self.output(x)

        return x * mask


class Predictor(nn.Module):
    def __init__(self, channels=64):
        super().__init__()

        self.model = CNN(channels)

        filter_kernel = gaussian_kernel(kernel_size=7, std=1.0)

        self.register_buffer(name="filter_kernel", tensor=filter_kernel)

    def forward(self, x, mask):
        B, N, C = x.shape

        mask = mask.view(B, 1, N, 1)

        x = x.transpose(1, 2).view(B, C, N, 1)

        x = self.model(x, mask)

        x = x.view(B, 5, N)

        x = F.pad(x, pad=(3, 3), mode="constant", value=0.0)

        x = x.unfold(dimension=2, size=7, step=1)

        x = torch.einsum("bcnm,m->bcn", x, self.filter_kernel)

        return x


def seed_all(seed):
    random.seed(seed)
    numpy.random.seed(seed)
    torch.manual_seed(seed)

    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True


def gaussian(x, std):
    pi = torch.tensor(math.pi)
    s2 = 2.0 * torch.tensor(std).square()
    x2 = torch.tensor(x).square().neg()

    return torch.exp(x2 / s2) * torch.rsqrt(s2 * pi)


def gaussian_kernel(kernel_size, std=1.0):
    kernel = [gaussian(i - (kernel_size // 2), std) for i in range(kernel_size)]

    kernel = torch.tensor(kernel)
    kernel = kernel / kernel.sum()

    return kernel


def make_mask(embeddings):
    B, N, _ = embeddings.shape

    mask = torch.zeros((B, N), dtype=embeddings.dtype, device=embeddings.device)
    mask[0, :N] = 1.0

    return mask


def load_models(device):
    model_files = [
        "cv_0.pt",
        "cv_1.pt",
        "cv_2.pt",
        "cv_3.pt",
        "cv_4.pt",
    ]
    models = []

    for model_file in model_files:
        model = Predictor()
        file_path = os.path.join(
            os.path.dirname(__file__), "data", "tmbed_models", model_file
        )
        log.debug(f"Loading model {file_path}")
        model.load_state_dict(torch.load(file_path)["model"])

        model = model.eval().to(device)

        models.append(model)

    return models


def predict_sequences(models, embedding, mask):
    B, N, _ = embedding.shape

    num_models = len(models)

    with torch.no_grad():
        pred = torch.zeros((B, 5, N), device=embedding.device)

        for model in models:
            y = model(embedding, mask)
            pred = pred + torch.softmax(y, dim=1)

        pred = pred / num_models

    return pred.detach()


def tmbed(sequences, models_path):
    log.info("Predicting transmembrane domains...")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.debug(f"Device is {device}")

    log.debug("Loading encoder")
    encoder = T5Encoder(models_path, torch.cuda.is_available())
    log.debug("Loading decoder")
    decoder = Decoder()

    models = load_models(device)

    pred_map = {0: "B", 1: "b", 2: "H", 3: "h", 4: "S", 5: "i", 6: "o"}

    states = {
        "B": "beta_outwards",
        "b": "beta_inwards",
        "H": "alpha_outwards",
        "h": "alpha_inwards",
        "S": "signal_peptide",
        "i": "inside",
        "o": "outside",
    }

    for sequence in sequences:
        log.debug
        try:
            log.debug(f"Predicting transmembrane domains for {sequence.id}...")
            embedding = encoder.embed(sequence.seq)
        except torch.cuda.OutOfMemoryError:
            log.warning(
                f"GPU ran out of memory when encoding {sequence.id} - skipping..."
            )
            continue

        # CPU alternative, implement fallback?
        # encoder.to_cpu()
        # torch.cuda.empty_cache()
        # embeddings = encoder.embed(sequences)

        embedding = embedding.to(device=device)
        embedding = embedding.to(dtype=torch.float32)

        mask = make_mask(embedding)

        probabilities = predict_sequences(models, embedding, mask)

        mask = mask.cpu()
        probabilities = probabilities.cpu()

        prediction = decoder(probabilities, mask).byte()

        sequence.transmembrane_predictions = "".join(
            pred_map[int(x)] for x in prediction[0, : len(sequence.seq)]
        )

        # Annotate relevant transmembrane domains
        for i, state in enumerate(sequence.transmembrane_predictions):
            if i == 0:
                previous_state = state
                state_start = i
                continue

            if state != previous_state:
                sequence.add_annotation(
                    states[previous_state],
                    "tmbed",
                    state_start + 1,
                    i,
                )

                state_start = i
                previous_state = state

        # Add final annotation
        sequence.add_annotation(
            states[previous_state],
            "tmbed",
            state_start + 1,
            len(sequence.seq),
        )

    return sequences
