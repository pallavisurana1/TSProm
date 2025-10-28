
# Usage instructions:
# conda activate dnabert2.0
# # cd /data/private/psurana/ramanaServer_code/scripts/OtherModel_FineTune/final
# # CUDA_VISIBLE_DEVICES=0,1,2 TOKENIZERS_PARALLELISM=false PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True nohup python aug_2025_genalm_1.py > finetune_all.log 2>&1 &

# cd /data/private/psurana/TSProm/src/1_finetune/
# CUDA_VISIBLE_DEVICES=5,6,7 TOKENIZERS_PARALLELISM=false PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True python FineTune_GENALM.py

# export NEPTUNE_PROJECT="lavi.ana/GenaLM-human-mouse-Finetune"
# export TOKENIZERS_PARALLELISM=false
# export NEPTUNE_API_TOKEN="eyJhcGlfYWRkcmVzcyI6Imh0dHBzOi8vYXBwLm5lcHR1bmUuYWkiLCJhcGlfdXJsIjoiaHR0cHM6Ly9hcHAubmVwdHVuZS5haSIsImFwaV9rZXkiOiIwNzg0NzY4MC1lNzQ5LTRkZDUtYTAwYi0zNDE4YzZjNzk0MTQifQ=="


def run_genalm_finetune(data_base_path, output_base_path, tissue_folder, length,
                        learning_rate=2e-5, batch_size=12, max_length=2000, num_train_epochs=10):
    """
    Fine-tune GENA-LM on tissue-specific promoter data.

    Parameters:
    - data_base_path (str): Root directory containing all data folders
    - output_base_path (str): Where to save model and evaluation results
    - tissue_folder (str): e.g., "tsp_brain_low"
    - length (str): e.g., "2000", "3000", "4000"
    - learning_rate (float): Training learning rate
    - batch_size (int): Batch size per device
    - max_length (int): Max token length for tokenizer
    - num_train_epochs (int): Number of training epochs
    """
    
    import torch
    import os
    import pandas as pd
    import json
    from transformers import AutoTokenizer, AutoModel, TrainingArguments, Trainer
    from transformers import DataCollatorWithPadding
    from datasets import Dataset, DatasetDict
    import numpy as np
    from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, matthews_corrcoef
    import importlib
    
    # Debug print
    print(f"🔧 Running with: LR={learning_rate}, Batch={batch_size}, Epochs={num_train_epochs}, Max_len={max_length}")
    print(f"📁 Tissue: {tissue_folder} | Length: {length}")
    print(f"🔢 CUDA devices available: {torch.cuda.device_count()}")
    
    # Model selection based on sequence length
    if length == "2000":
        tokenizer = AutoTokenizer.from_pretrained('AIRI-Institute/gena-lm-bert-large-t2t')
        model_name = 'AIRI-Institute/gena-lm-bert-large-t2t'
        model_class_name = "BertForSequenceClassification"
    else:  # For 3000 and 4000
        tokenizer = AutoTokenizer.from_pretrained('AIRI-Institute/gena-lm-bigbird-base-t2t')
        model_name = 'AIRI-Institute/gena-lm-bigbird-base-t2t'
        model_class_name = "BigBirdForSequenceClassification"

    # Fix padding token if needed
    if tokenizer.pad_token is None:
        tokenizer.pad_token = tokenizer.eos_token

    # Load base model to get module info
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
    gena_module_name = model.__class__.__module__
    print(f"📦 Using module: {gena_module_name}")

    # Get classification model class
    cls = getattr(importlib.import_module(gena_module_name), model_class_name)
    
    # Create classification model
    model = cls.from_pretrained(model_name, num_labels=2)
    print(f'🧠 Classification head: {model.classifier}')

    # Setup paths
    trial_data = os.path.join(data_base_path, length, tissue_folder)
    output_finish_path = os.path.join(output_base_path, tissue_folder, length, f"lr_{learning_rate:.0e}")
    
    # Create output directory
    os.makedirs(output_finish_path, exist_ok=True)
    print(f"📂 Output directory: {output_finish_path}")

    # Load data with debugging
    try:
        dev = pd.read_csv(os.path.join(trial_data, 'dev.csv'))
        train = pd.read_csv(os.path.join(trial_data, 'train.csv'))
        test = pd.read_csv(os.path.join(trial_data, 'test.csv'))
        
        # Debug: Check column names and data types
        print(f"📊 Train columns: {list(train.columns)}")
        print(f"🔍 Sample train data:\n{train.head(2)}")
        print(f"🔍 Sequence column type: {type(train['Sequence'].iloc[0])}")
        
        print(f"📊 Data loaded - Train: {len(train)}, Dev: {len(dev)}, Test: {len(test)}")
    except Exception as e:
        print(f"❌ Error loading data: {e}")
        raise

    # Clean and prepare data - remove pandas index column that might cause issues
    train = train.reset_index(drop=True)
    dev = dev.reset_index(drop=True)
    test = test.reset_index(drop=True)
    
    # Convert to the format expected by the model
    def prepare_data(df):
        return {
            'sequence': df['Sequence'].astype(str).tolist(),
            'label': df['Label'].astype(int).tolist()
        }
    
    train_data = prepare_data(train)
    dev_data = prepare_data(dev)
    test_data = prepare_data(test)
    
    # Create datasets from dictionaries
    train_dataset = Dataset.from_dict(train_data)
    dev_dataset = Dataset.from_dict(dev_data)
    test_dataset = Dataset.from_dict(test_data)

    dataset = DatasetDict({
        "train": train_dataset,
        "validation": dev_dataset,
        "test": test_dataset,
    })
    
    # Debug: Check first example
    print("🔍 First train example:", dataset["train"][0])

    # Tokenization - simpler approach
    def preprocess_function(examples):
        # Tokenize the sequences
        tokenized = tokenizer(
            examples["sequence"],
            truncation=True, 
            max_length=max_length,
            padding=False,  # Let DataCollator handle padding
            return_tensors=None
        )
        # Keep the labels
        tokenized["labels"] = examples["label"]  # Use "labels" not "label"
        return tokenized

    tokenized_dataset = dataset.map(preprocess_function, batched=True)
    
    # Remove original columns to avoid conflicts
    tokenized_dataset = tokenized_dataset.remove_columns(["sequence", "label"])
    
    # Debug: Check tokenized data
    print("🔍 Tokenized sample:", {k: (type(v), len(v) if isinstance(v, list) else v) for k, v in tokenized_dataset["train"][0].items()})
    # Metrics computation
    def compute_metrics(eval_pred):
        predictions, labels = eval_pred
        predictions = np.argmax(predictions, axis=1)
        
        return {
            "accuracy": accuracy_score(labels, predictions),
            "precision": precision_score(labels, predictions, zero_division=0),
            "recall": recall_score(labels, predictions, zero_division=0),
            "f1": f1_score(labels, predictions, zero_division=0),
            "matthews_correlation": matthews_corrcoef(labels, predictions)
        }

    # Create data collator
    data_collator = DataCollatorWithPadding(tokenizer=tokenizer, return_tensors="pt")

    # Training arguments
    training_args = TrainingArguments(
        output_dir=output_finish_path,
        learning_rate=learning_rate,
        lr_scheduler_type="constant_with_warmup",
        warmup_ratio=0.1,
        optim='adamw_torch',
        weight_decay=0.0,
        per_device_train_batch_size=batch_size,
        per_device_eval_batch_size=batch_size,
        num_train_epochs=num_train_epochs,
        evaluation_strategy="epoch",
        save_strategy="epoch",
        logging_strategy="epoch",
        load_best_model_at_end=True,
        dataloader_pin_memory=False,
        remove_unused_columns=False,
        report_to=None,
        # Added for better multi-GPU support
        dataloader_num_workers=4,
        fp16=True,  # Enable mixed precision for better memory usage
    )

    # Create trainer
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=tokenized_dataset["train"],
        eval_dataset=tokenized_dataset["validation"],
        tokenizer=tokenizer,
        data_collator=data_collator,
        compute_metrics=compute_metrics,
    )

    print("🚀 Starting training...")
    trainer.train()

    # Evaluate on test set
    print("📈 Evaluating on test set...")
    results = trainer.evaluate(eval_dataset=tokenized_dataset["test"])
    print("📊 Test Results:")
    for key, value in results.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.4f}")
        else:
            print(f"  {key}: {value}")

    # Save results
    output_json_path = os.path.join(output_finish_path, "eval_results.json")
    with open(output_json_path, "w") as f:
        json.dump(results, f, indent=4)
    print(f"💾 Results saved to: {output_json_path}")

    # Cleanup memory
    del model, trainer
    torch.cuda.empty_cache()
    print("Memory cleaned up")


# ============================================================================
# Main execution loop with error handling
# ============================================================================

# tissue_dirs = [
#     "tsp_brain_low", "tsp_brain_tenh", "tsp_liver_low", "tsp_liver_tenh",
#     "tsp_testis_low", "tsp_brain_null", "tsp_brain_wide",
#     "tsp_liver_null", "tsp_liver_wide", "tsp_testis_null", "tsp_testis_wide", "tsp_testis_tenh"
# ]
# tissue_dirs =["tsp_spleen_null", "tsp_spleen_wide", "tsp_spleen_tenh", "tsp_spleen_low", 
#               "tsp_muscle_null", "tsp_muscle_wide"]

tissue_dirs = ["tsp_muscle_tenh", "tsp_muscle_low"]
lengths = ["3000", "4000", "2000"]
length_map = {"2000": 2000, "3000": 3000, "4000": 4000}
learning_rates = [2e-5, 3e-4, 5e-6]

import os

print("🎯 Starting GENA-LM fine-tuning experiments...")
print(f"📋 Total experiments: {len(tissue_dirs)} × {len(lengths)} × {len(learning_rates)} = {len(tissue_dirs) * len(lengths) * len(learning_rates)}")

experiment_count = 0
total_experiments = len(tissue_dirs) * len(lengths) * len(learning_rates)

for tissue in tissue_dirs:
    for length in lengths:
        for lr in learning_rates:
            experiment_count += 1
            batch_size = 8 if "testis" in tissue else 12
            lr_tag = f"lr_{lr:.0e}"
            output_dir = os.path.join(
                "/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/gena-lm_finetune/unique_TSS/human_mouse",
                tissue, length, lr_tag
            )
            
            print(f"\n{'='*80}")
            print(f"🧪 Experiment {experiment_count}/{total_experiments}")
            print(f"🧬 Tissue: {tissue} | Length: {length} | LR: {lr} | Batch: {batch_size}")
            print(f"{'='*80}")
            
            if os.path.exists(os.path.join(output_dir, "eval_results.json")):
                print(f"✅ Skipping: {tissue} | Length: {length} | LR: {lr} (Already completed)")
                continue

            try:
                run_genalm_finetune(
                    data_base_path="/data/projects/dna/pallavi/data_TSp_Vs_Rest",
                    output_base_path="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/gena-lm_finetune/unique_TSS/human_mouse",
                    tissue_folder=tissue,
                    length=length,
                    learning_rate=lr,
                    batch_size=batch_size,
                    max_length=length_map[length],
                    num_train_epochs=10
                )
                print(f"✅ Completed: {tissue} | Length: {length} | LR: {lr}")
                
            except Exception as e:
                print(f"❌ Error for {tissue} | Length: {length} | LR: {lr}: {str(e)}")
                import traceback
                print(f"📝 Full traceback:\n{traceback.format_exc()}")
                continue

print(f"\n🎉 All experiments completed! Check logs for any errors.")