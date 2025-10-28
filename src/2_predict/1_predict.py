#!/usr/bin/env python3
import os, csv, torch, numpy as np
import transformers
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Sequence
from torch.utils.data import Dataset, DataLoader, SequentialSampler
import torch.nn.functional as F
from tqdm import tqdm

# ──────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────
@dataclass
class ModelArguments:
    model_name_or_path: Optional[str] = field(default=None)


@dataclass
class DataArguments:
    data_path: str = field(default=None, metadata={"help": "Directory with train.csv/dev.csv/test.csv"})
    kmer: int = field(default=-1)


@dataclass
class PredictionArguments(transformers.TrainingArguments):
    cache_dir: Optional[str] = field(default=None)
    output_dir: str = field(default="pred_output")
    per_device_eval_batch_size: int = field(default=8)
    model_max_length: int = field(default=512)
    save_all: bool = field(default=False, metadata={"help": "Save all predictions, not just correct ones."})


# ──────────────────────────────────────────────
# Dataset
# ──────────────────────────────────────────────
class SupervisedDataset(Dataset):
    def __init__(self, data_path: str, tokenizer: transformers.PreTrainedTokenizer, kmer: int = -1):
        super().__init__()

        with open(data_path, "r") as f:
            data = list(csv.reader(f))[1:]

        if len(data[0]) == 2:
            texts = [d[0] for d in data]
            labels = [int(d[1]) for d in data]
        else:
            raise ValueError(f"{data_path} must have [sequence,label] columns")

        if kmer != -1:
            def generate_kmer_str(sequence, k):
                return " ".join([sequence[i:i+k] for i in range(len(sequence) - k + 1)])
            texts = [generate_kmer_str(seq, kmer) for seq in texts]

        enc = tokenizer(
            texts,
            return_tensors="pt",
            padding="longest",
            max_length=tokenizer.model_max_length,
            truncation=True,
        )
        self.input_ids = enc["input_ids"]
        self.attention_mask = enc["attention_mask"]
        self.labels = labels
        self.texts = texts  # Save original texts for decoding

    def __len__(self): return len(self.input_ids)
    def __getitem__(self, i) -> Dict[str, torch.Tensor]:
        return dict(
            input_ids=self.input_ids[i],
            attention_mask=self.attention_mask[i],
            labels=torch.tensor(self.labels[i])
        )


class DataCollatorForSupervisedDataset:
    def __init__(self, tokenizer): self.tokenizer = tokenizer
    def __call__(self, batch: Sequence[Dict]) -> Dict[str, torch.Tensor]:
        input_ids = [b["input_ids"] for b in batch]
        labels = [b["labels"] for b in batch]

        input_ids = torch.nn.utils.rnn.pad_sequence(
            input_ids, batch_first=True, padding_value=self.tokenizer.pad_token_id
        )
        attn_mask = input_ids.ne(self.tokenizer.pad_token_id)

        return dict(
            input_ids=input_ids,
            attention_mask=attn_mask,
            labels=torch.tensor(labels).long()
        )


# ──────────────────────────────────────────────
# Prediction function
# ──────────────────────────────────────────────
def run_prediction_for_split(split, file_path, model, tokenizer, pred_args, kmer, device):
    dataset = SupervisedDataset(file_path, tokenizer, kmer)
    collator = DataCollatorForSupervisedDataset(tokenizer)

    sampler = SequentialSampler(dataset)
    dataloader = DataLoader(dataset, sampler=sampler, batch_size=pred_args.per_device_eval_batch_size, collate_fn=collator)

    if torch.cuda.device_count() > 1 and not isinstance(model, torch.nn.DataParallel):
        model = torch.nn.DataParallel(model)

    model.eval()
    preds, true_labels = [], []

    for batch in tqdm(dataloader, desc=f"Predicting {split}"):
        input_ids = batch["input_ids"].to(device)
        attn_mask = batch["attention_mask"].to(device)
        labels = batch["labels"].to(device)

        with torch.no_grad():
            outputs = model(input_ids=input_ids, attention_mask=attn_mask, labels=labels)
            logits = outputs.logits

        preds.append(logits.detach().cpu())
        true_labels.append(labels.detach().cpu())

    preds = torch.cat(preds, dim=0)
    true_labels = torch.cat(true_labels, dim=0)
    probs = F.softmax(preds, dim=1).numpy()
    pred_labels = np.argmax(probs, axis=1)
    true_labels = true_labels.numpy()

    all_csv = os.path.join(pred_args.output_dir, f"{split}_all_predictions.csv")
    os.makedirs(pred_args.output_dir, exist_ok=True)

    with open(all_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sequence", "true_label", "pred_label", "prob_0", "prob_1"])
        for i, seq in enumerate(dataset.input_ids):
            writer.writerow([
                tokenizer.decode(seq, skip_special_tokens=True),
                int(true_labels[i]),
                int(pred_labels[i]),
                probs[i][0],
                probs[i][1]
            ])

    print(f"✅ All predictions for {split} saved to {all_csv}")


def main():
    parser = transformers.HfArgumentParser((ModelArguments, DataArguments, PredictionArguments))
    model_args, data_args, pred_args = parser.parse_args_into_dataclasses()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    tokenizer = transformers.AutoTokenizer.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=pred_args.cache_dir,
        model_max_length=pred_args.model_max_length,
        use_fast=True,
        trust_remote_code=True,
    )
    model = transformers.AutoModelForSequenceClassification.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=pred_args.cache_dir,
        trust_remote_code=True,
    ).to(device)

    for split in ["train", "dev", "test"]:
        file_path = os.path.join(data_args.data_path, f"{split}.csv")
        if os.path.exists(file_path):
            run_prediction_for_split(split, file_path, model, tokenizer, pred_args, data_args.kmer, device)


if __name__ == "__main__":
    main()
