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
import torch

from tqdm import tqdm
from pathlib import Path

from .model import Predictor
from .embed import T5Encoder
from .viterbi import Decoder

from .utils import seed_all, make_batches, collate_batch, make_mask


def load_models(device):
    models = []
    root_path = Path(__file__).parent
    model_path = Path(root_path, 'models/cnn/')

    for model_file in sorted(model_path.glob('*.pt')):
        model = Predictor()

        model.load_state_dict(torch.load(model_file)['model'])

        model = model.eval().to(device)

        models.append(model)

    return models


def predict_sequences(models, embeddings, mask):
    B, N, _ = embeddings.shape

    num_models = len(models)

    with torch.no_grad():
        pred = torch.zeros((B, 5, N), device=embeddings.device)

        for model in models:
            y = model(embeddings, mask)
            pred = pred + torch.softmax(y, dim=1)

        pred = pred / num_models

    return pred.detach()

def tmbed(sequences):
    batch_size = 4000

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    encoder = T5Encoder("resistify_models/prott5/", torch.cuda.is_available())
    decoder = Decoder()

    models = load_models(device)

    predictions = dict()

    # Need to sort sequences by length? Why?
    sequences = sorted(sequences, key=lambda seq: len(seq.length))

    batches = make_batches(sequences, batch_size)

    for a, b in batches:
        batch = sequences[a:b]

        lengths = [len(sequence.sequence) for sequence in batch]
        sequences = [sequence.sequence for sequence in batch]

        embeddings = encoder.embed(sequences)
        # CPU alternative, implement fallback?
        #encoder.to_cpu()
        #torch.cuda.empty_cache()
        #embeddings = encoder.embed(sequences)

        embeddings = embeddings.to(device=device)
        embeddings = embeddings.to(dtype=torch.float32)

        mask = make_mask(embeddings, lengths)

        probabilities = predict_sequences(models, embeddings, mask)

        # Why?
        mask = mask.cpu()
        probabilities = probabilities.cpu()

        prediction = decoder(probabilities, mask).byte()

        for idx, sequence in enumerate(batch):
            length = len(sequence.sequence)
            for pred in prediction[idx, :length]:
                print(pred)
