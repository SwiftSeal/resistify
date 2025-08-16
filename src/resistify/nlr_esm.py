import torch
from transformers import AutoModelForTokenClassification
from rich.progress import Progress
from resistify.annotation import Protein, Annotation
from resistify.console import console
import logging

# Hide the stupid weights error
logging.getLogger("transformers.modeling_utils").setLevel(logging.ERROR)

def predict_lrr(proteins: list[Protein], device: str, batch_size: int = 32, lrr_gap: int = 75, lrr_length: int = 4):
    with Progress(console=console) as progress:
        task = progress.add_task("Predicting LRR motifs", total=len(proteins))

        model_name = "MoraySmith/LRR_ESMplusplus_small"
        console.log("Loading LRR prediction model")
        model = AutoModelForTokenClassification.from_pretrained(model_name, trust_remote_code=True).to(device)
        tokenizer = model.tokenizer

        for i in range(0, len(proteins), batch_size):
            batch_proteins = proteins[i:i + batch_size]
            batch_sequences = [p.sequence for p in batch_proteins]

            inputs = tokenizer(batch_sequences, return_tensors="pt", padding=True, truncation=True).to(device)
            with torch.no_grad():
                outputs = model(**inputs)

            predictions = torch.argmax(outputs.logits, dim=-1)

            # Iterate through the predictions and proteins in the current batch
            for protein, prediction_tensor, attention_mask_tensor in zip(batch_proteins, predictions, inputs['attention_mask']):
                labels = prediction_tensor.tolist()
                attention_mask = attention_mask_tensor.tolist()
                # Remove padding tokens
                labels = [label for label, mask in zip(labels, attention_mask) if mask == 1]
                labels = labels[1:-1]  # Remove CLS and SEP tokens
                if protein.length != len(labels):
                    raise ValueError(
                        f"Sequence length {len(protein.sequence)} does not match label length {len(labels)} for protein {protein.id}."
                    )

                motif_positions = [j for j, label in enumerate(labels) if label == 1]

                if motif_positions:
                    start, end, count = 0, 0, 0
                    for pos in motif_positions:
                        # Add LRR motif (make position 1-indexed)
                        protein.add_annotation(
                            Annotation(
                                "LxxLxL",
                                pos + 1,
                                pos + 6,
                                type="motif"
                            )
                        )
                        if pos - end < lrr_gap:
                            # Extend if still within gap limit
                            end = pos
                            count += 1
                        else:
                            if count >= lrr_length:
                                # Annotate if big enough
                                protein.add_annotation(
                                    Annotation(
                                        "LRR",
                                        start + 1,
                                        end + 1,
                                    )
                                )
                            # Then reset the counters
                            start = pos
                            end = pos
                            count = 0
                    if count >= lrr_length:
                        # Edge case for the last motif
                        protein.add_annotation(
                            Annotation(
                                "LRR",
                                start + 1,
                                end + 1,
                            )
                        )
                progress.update(task, advance=1)