#! /usr/bin/env python3
"""
This script is a simple example of how to fine-tune a Synthyra FastPLM model for a protein sequence regression or classification task.
For regression we look at the binding affinity of two proteins (pkd)
For classification we look at the solubility of a protein (membrane bound or not)
"""
import numpy as np
from sklearn.model_selection import train_test_split
from datasets import Dataset
from typing import Any
from transformers import (
    AutoModelForTokenClassification,
    Trainer,
    TrainingArguments,
    EarlyStoppingCallback,
    DataCollatorForTokenClassification
)
from peft import LoraConfig, get_peft_model


# Shared arguments for the trainer
BASE_TRAINER_KWARGS = {
    "warmup_steps": 500,
    "weight_decay": 0.01,
    "logging_steps": 100,
    "eval_strategy": "steps",
    "eval_steps": 500,
    "save_strategy": "steps",
    "save_steps": 500,
    "load_best_model_at_end": True,
    "metric_for_best_model": "eval_loss",
    "greater_is_better": False,
    "report_to": "none",
    "label_names": ["labels"]
}

def prepare_lrr_data(tokenizer):
    sequences = []
    labels = []

    with open("training/data/LRRPredictor_validation_groundTruth.csv") as f:
        for line in f:
            _, data_type, value = line.strip().split(",")
            if data_type == "sequence":
                sequences.append(value)
            elif data_type == "LRR motifs":
                label = np.zeros(len(value), dtype=int)
                for i, c in enumerate(value):
                    if c == "-":
                        label[i] = 0
                    else:
                        label[i] = 1
                labels.append(label)

    train_sequences, test_sequences, train_labels, test_labels = train_test_split(sequences, labels, test_size=0.25, shuffle=True)

    train_tokenized = tokenizer(train_sequences)
    test_tokenized = tokenizer(test_sequences)

    train_dataset = Dataset.from_dict(train_tokenized)
    test_dataset = Dataset.from_dict(test_tokenized)

    train_dataset = train_dataset.add_column("labels", train_labels)
    test_dataset = test_dataset.add_column("labels", test_labels)

    return train_dataset, test_dataset

def prepare_cc_data(tokenizer):
    sequences = []
    labels = []

    with open("training/data/deepcoil_train.csv") as f:
        next(f)  # Skip header
        for line in f:
            _, seq, cc = line.strip().split(",")
            assert len(seq) == len(cc), f"Length mismatch for {id}: {len(seq)} vs {len(cc)}"
            sequences.append(seq)
            label = np.zeros(len(seq), dtype=int)
            for i, c in enumerate(cc):
                if c == '1':
                    label[i] = 1
                elif c == '0':
                    label[i] = 0
                else:
                    raise ValueError(f"Unexpected label character: {c}")
            labels.append(label)
            assert len(label) == len(seq), f"Label length mismatch for {id}: {len(label)} vs {len(seq)}"
    
    train_sequences, test_sequences, train_labels, test_labels = train_test_split(sequences, labels, test_size=0.25, shuffle=True)

    train_tokenized = tokenizer(train_sequences)
    test_tokenized = tokenizer(test_sequences)

    train_dataset = Dataset.from_dict(train_tokenized)
    test_dataset = Dataset.from_dict(test_tokenized)

    train_dataset = train_dataset.add_column("labels", train_labels)
    test_dataset = test_dataset.add_column("labels", test_labels)

    return train_dataset, test_dataset

# Get the model ready, with or without LoRA
def initialize_model(model_name: str, num_labels: int, use_lora: bool = True, lora_config: Any = None):
    """
    Initialize a model with optional LoRA support
    
    Args:
        model_name: Name or path of the pretrained model
        num_labels: Number of labels for the task (1 for regression)
        use_lora: Whether to use LoRA for fine-tuning
        lora_config: Custom LoRA configuration (optional)
        
    Returns:
        model: The initialized model
        tokenizer: The model's tokenizer
    """
    print(f"Loading model {model_name} with {num_labels} labels...")
    
    # Load base model
    model = AutoModelForTokenClassification.from_pretrained(
        model_name,
        trust_remote_code=True,
        num_labels=num_labels
    )
    tokenizer = model.tokenizer
    
    # Apply LoRA if requested
    if use_lora:
        # Default LoRA configuration if none provided
        if lora_config is None:
            # Target modules for ESM++ or ESM2 models
            target_modules = ["layernorm_qkv.1", "out_proj", "query", "key", "value", "dense"]
            
            lora_config = LoraConfig(
                r=8,
                lora_alpha=16,
                lora_dropout=0.01,
                bias="none",
                target_modules=target_modules,
            )
        
        # Apply LoRA to the model
        model = get_peft_model(model, lora_config)
        
        # Unfreeze the classifier head
        for param in model.classifier.parameters():
            param.requires_grad = True
        
        # Print parameter statistics
        total_params = sum(p.numel() for p in model.parameters())
        trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
        non_trainable_params = total_params - trainable_params
        print(f"Total parameters: {total_params}")
        print(f"Trainable parameters: {trainable_params}")
        print(f"Non-trainable parameters: {non_trainable_params}")
        print(f"Percentage of parameters being trained: {100 * trainable_params / total_params:.2f}%")
    
    return model, tokenizer


# Main function
if __name__ == "__main__":
    print("Loading datasets for classification task...")

    model_name = 'Synthyra/ESMplusplus_small'
    num_labels = 2
    use_lora = False
    num_epochs = 10
    gradient_accumulation_steps = 1
    batch_size = 8
    learning_rate = 5e-5
    patience = 3

    # Initialize model with modular function
    model, tokenizer = initialize_model(
        model_name=model_name,
        num_labels=num_labels,
        use_lora=use_lora,
    )

    train_dataset, test_dataset = prepare_lrr_data(tokenizer)

    # Create data collator
    data_collator = DataCollatorForTokenClassification(tokenizer=tokenizer)
    
    # Define training arguments
    output_dir = "./results_classification_lora" if use_lora else "./results_classification"
    logging_dir = "./logs_classification_lora" if use_lora else "./logs_classification"
    
    training_args = TrainingArguments(
        output_dir=output_dir,
        num_train_epochs=num_epochs,
        gradient_accumulation_steps=gradient_accumulation_steps,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        logging_dir=logging_dir,
        learning_rate=learning_rate,
        **BASE_TRAINER_KWARGS
    )
    
    # Create trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        eval_dataset=test_dataset,
        data_collator=data_collator,
        #callbacks=[EarlyStoppingCallback(early_stopping_patience=patience)]
    )

    trainer.train()
