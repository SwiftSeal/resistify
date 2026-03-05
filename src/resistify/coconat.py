import torch
import esm
import tqdm.auto as tqdm
from transformers import T5EncoderModel, T5Tokenizer
import torch.nn as nn
import logging
from pathlib import Path
from resistify.annotation import Protein

logger = logging.getLogger(__name__)


class TransposeX(nn.Module):
    def __init__(self):
        super(TransposeX, self).__init__()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
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

    def forward(self, x: torch.Tensor, lengths: list[int]) -> torch.Tensor:
        x = self.cnn(x)

        len_tensor = torch.as_tensor(lengths, dtype=torch.int64, device="cpu")

        x = torch.nn.utils.rnn.pack_padded_sequence(
            x, len_tensor, batch_first=True, enforce_sorted=False
        )
        x, _ = self.lstm(x)
        x, _ = torch.nn.utils.rnn.pad_packed_sequence(x, batch_first=True)

        x = self.linear(x)
        return self.final(x)


class CoCoNatPredictor:
    def __init__(self, device: str):
        self.device = torch.device(device)
        self._load_models()

    def _load_models(self):
        self.t5_tokenizer = T5Tokenizer.from_pretrained(
            "Rostlab/prot_t5_xl_half_uniref50-enc", do_lower_case=False, legacy=True
        )
        self.t5_model = (
            T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
            .to(self.device)
            .eval()
        )

        self.esm_model, self.esm_alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.esm_model = self.esm_model.to(self.device).eval()
        self.esm_batch_converter = self.esm_alphabet.get_batch_converter()

        model_path = Path(__file__).parent / "data" / "coconat.ckpt"
        checkpoint = torch.load(model_path, map_location=self.device)
        self.reg_model = MMModelLSTM().to(self.device).eval()
        self.reg_model.load_state_dict(checkpoint["state_dict"])

    @torch.inference_mode()
    def predict_batch(self, batch_proteins: list[tuple[str, str]]):
        ids, seqs = zip(*batch_proteins)
        lengths = [len(s) for s in seqs]

        t5_seqs = [
            " ".join(
                list(
                    s.replace("U", "X")
                    .replace("Z", "X")
                    .replace("O", "X")
                    .replace("B", "X")
                )
            )
            for s in seqs
        ]
        t5_inputs = self.t5_tokenizer(t5_seqs, return_tensors="pt", padding=True).to(
            self.device
        )
        t5_emb = self.t5_model(**t5_inputs).last_hidden_state

        t5_res = [t5_emb[i, : lengths[i]] for i in range(len(lengths))]

        _, _, esm_tokens = self.esm_batch_converter(batch_proteins)
        esm_out = self.esm_model(esm_tokens.to(self.device), repr_layers=[33])[
            "representations"
        ][33]
        esm_res = [esm_out[i, 1 : lengths[i] + 1] for i in range(len(lengths))]

        results = {}
        for i in range(len(lengths)):
            merged = torch.cat([t5_res[i], esm_res[i]], dim=-1).unsqueeze(0)
            pred = self.reg_model(merged, [lengths[i]])
            cc_prob = 1 - pred[0, :, 0]
            results[ids[i]] = cc_prob.cpu().numpy().tolist()

        return results


def predict_coils(
    proteins: dict[str, Protein], device: str, batch_size: int, threads: int
):
    logger.info("Running CoCoNat to predict coiled-coils")
    torch.set_num_threads(threads)
    predictor = CoCoNatPredictor(device=device)

    to_process = []
    for p in proteins.values():
        if p.nbarc_start and p.nbarc_start > 5:
            seq = p.sequence[: p.nbarc_start - 1]
            to_process.append((p.id, seq))

    batches = [
        to_process[i : i + batch_size] for i in range(0, len(to_process), batch_size)
    ]
    for batch in tqdm.tqdm(batches, desc="CoCoNat", unit="batch"):
        batch_results = predictor.predict_batch(batch)

        for p_id, cc_probs in batch_results.items():
            proteins[p_id].cc_probs = cc_probs
            proteins[p_id].annotate_cc()
    logger.info("CoCoNat predictions completed")
