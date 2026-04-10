import os
import math
import torch
import torch.nn as nn
import torch.nn.functional as F
import logging
from tqdm.auto import tqdm
from resistify.annotation import Protein, Annotation
from resistify.esm import ESM2Encoder, ESM_DIM

logger = logging.getLogger(__name__)

PRED_MAP = {0: "H", 1: "h", 2: "S", 3: "i", 4: "o"}

STATES = {
    "H": "alpha_outwards",
    "h": "alpha_inwards",
    "S": "signal_peptide",
    "i": "inside",
    "o": "outside",
}


class Decoder(nn.Module):
    def __init__(self):
        super().__init__()

        self._init_transitions()

    def _init_transitions(self):
        num_tags = 17

        end_transitions = torch.full((num_tags,), -100)
        start_transitions = torch.full((num_tags,), -100)

        transitions = torch.full((num_tags, num_tags), -100)

        for i in [0, 5, 10, 15, 16]:
            start_transitions[i] = 0  # H1a, H1b, S1, i, o

        for i in range(4):
            transitions[0 + i, 1 + i] = 0   # Hxa -> Hya
            transitions[5 + i, 6 + i] = 0   # Hxb -> Hyb
            transitions[10 + i, 11 + i] = 0  # Sx -> Sy

        for i in [4, 9, 14]:
            transitions[i, i] = 0  # X5 -> X5

        transitions[4, 16] = 0    # H5a -> o
        transitions[9, 15] = 0    # H5b -> i
        transitions[14, 15:] = 0  # S5  -> (i, o)

        transitions[15, 0] = 0    # i -> H1a
        transitions[15, 15:] = 0  # i -> (i, o)

        transitions[16, 5] = 0    # o -> H1b
        transitions[16, 15:] = 0  # o -> (i, o)

        for i in [4, 9, 14, 15, 16]:
            end_transitions[i] = 0  # H5a, H5b, S5, i, o

        repeats = torch.tensor([10, 5, 1, 1], dtype=torch.int32)

        mapping = torch.arange(5, dtype=torch.int32)
        mapping = mapping.repeat_interleave(
            torch.tensor([5, 5, 5, 1, 1])  # H  # h  # S  # i  # o
        )

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

        nn.init.xavier_uniform_(self.conv.weight)

    def forward(self, x, mask):
        x = self.func(self.norm(self.conv(x), mask))

        return x * mask


class CNN(nn.Module):
    def __init__(self, channels):
        super().__init__()

        self.input = Conv(ESM_DIM, channels, 1, 1, 0)

        self.dwc1 = Conv(channels, channels, 9, 1, 4, groups=channels)
        self.dwc2 = Conv(channels, channels, 21, 1, 10, groups=channels)

        self.dropout = nn.Dropout2d(p=0.50, inplace=True)

        self.output = nn.Conv2d(3 * channels, 4, 1, 1, 0)

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

        x = x.view(B, 4, N)

        x = F.pad(x, pad=(3, 3), mode="constant", value=0.0)

        x = x.unfold(dimension=2, size=7, step=1)

        x = torch.einsum("bcnm,m->bcn", x, self.filter_kernel)

        return x


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


def make_mask(embeddings, attention_mask):
    return attention_mask.to(dtype=embeddings.dtype, device=embeddings.device)


def load_model(device):
    model = Predictor()
    file_path = os.path.join(
        os.path.dirname(__file__), "data", "tmbed_models", "tmbed_model.pt"
    )
    logger.debug(f"Loading model {file_path}")
    model.load_state_dict(torch.load(file_path, weights_only=True)["model"])
    return model.eval().to(device)


def predict_sequences(model, embedding, mask):
    B, N, _ = embedding.shape

    with torch.no_grad():
        pred = torch.softmax(model(embedding, mask), dim=1)

    return pred.detach()


def tmbed(proteins: dict[str, Protein], device: str, batch_size: int, threads: int, encoder: ESM2Encoder | None = None):
    torch.set_num_threads(threads)
    logger.info("Predicting transmembrane domains with TMBed")

    if encoder is None:
        encoder = ESM2Encoder(device)
    decoder = Decoder()
    model = load_model(device)

    protein_list = list(proteins.values())

    for i in tqdm(range(0, len(protein_list), batch_size)):
        batch = protein_list[i : i + batch_size]
        batch_seqs = [p.sequence for p in batch]

        try:
            embeddings, att_mask = encoder.embed(batch_seqs)
            embeddings = embeddings.to(device=device, dtype=torch.float32)

            mask = make_mask(embeddings, att_mask)
            probabilities = predict_sequences(model, embeddings, mask)

            predictions = decoder(probabilities.cpu(), mask.cpu()).byte()

            for idx, protein in enumerate(batch):
                seq_len = len(protein.sequence)
                pred_seq = predictions[idx, :seq_len]

                protein.transmembrane_predictions = "".join(
                    PRED_MAP[int(x)] for x in pred_seq
                )

                _annotate_protein(protein)

        except torch.cuda.OutOfMemoryError:
            logger.error("GPU OOM for batch - try reducing batch_size")
            raise


def _annotate_protein(protein):
    """Extracted annotation logic for readability."""
    for i, state in enumerate(protein.transmembrane_predictions):
        if i == 0:
            previous_state = state
            state_start = i
            continue

        if state != previous_state:
            if STATES[previous_state] not in ("inside", "outside"):
                protein.add_annotation(
                    Annotation(
                        name=STATES[previous_state],
                        start=state_start + 1,
                        end=i,
                        type="domain",
                        source="tmbed",
                    )
                )

            state_start = i
            previous_state = state

    # Add final annotation
    if STATES[previous_state] not in ("inside", "outside"):
        protein.add_annotation(
            Annotation(
                name=STATES[previous_state],
                start=state_start + 1,
                end=len(protein.sequence),
                type="domain",
                source="tmbed",
            )
        )
