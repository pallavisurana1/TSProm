#!/usr/bin/env python3
"""
DNABERT2 Attention Extraction + Motif Preprocessing (GPU-Flexible Version)
Date: 2025-08-29

Pipeline:
  1. Collect jobs (base_model, subset, lr_dir, checkpoint, data).
  2. Load DNABERT2 model + tokenizer with attention outputs.
  3. Extract attentions for train/dev/test splits.
  4. Expand into per-token records and compute base-pair bounds.
  5. Apply attention-based cutoff mask to highlight important tokens.
  6. Save results (Parquet + CSV for R compatibility).

Requires: torch, transformers, pandas, pyarrow, statsmodels, biopython, etc.

Run:
conda activate dnabert_aug_2025_jupyter
cd /data/private/psurana/DNA2-TSp/code/attention/; 
# export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:35000
# Run on specific GPU: CUDA_VISIBLE_DEVICES=2 python 1_raw_attention_extract-gpu.py
# Or let it auto-select: python 1_raw_attention_extract-gpu.py --gpu auto
# Or specify GPU ID: python 1_raw_attention_extract-gpu.py --gpu 1
python 1_raw_attention_extract-gpu-oneJob.py --gpu 0 --batch-size 256

"""

# ─────────────────────────────────────────────────────────────
# Imports
# ─────────────────────────────────────────────────────────────
import os, ast, re, glob, difflib, gc, argparse
import torch
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import ahocorasick
from itertools import combinations
from collections import Counter
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from Bio import pairwise2, motifs, SeqIO

# ─────────────────────────────────────────────────────────────
# Command Line Arguments
# ─────────────────────────────────────────────────────────────
def parse_args():
    parser = argparse.ArgumentParser(description='DNABERT2 Attention Extraction with GPU flexibility')
    parser.add_argument('--gpu', type=str, default='auto', 
                       help='GPU selection: "auto" (best available), "cpu", or specific GPU ID (e.g., "0", "1")')
    parser.add_argument('--batch-size', type=int, default=16,
                       help='Batch size for processing (adjust based on GPU memory)')
    parser.add_argument('--max-memory-fraction', type=float, default=0.8,
                       help='Maximum fraction of GPU memory to use')
    return parser.parse_args()

# ─────────────────────────────────────────────────────────────
# GPU Setup Functions
# ─────────────────────────────────────────────────────────────
def get_gpu_memory_info():
    """Get memory info for all available GPUs"""
    gpu_info = []
    if torch.cuda.is_available():
        for i in range(torch.cuda.device_count()):
            props = torch.cuda.get_device_properties(i)
            total_mem = props.total_memory / (1024**3)  # GB
            
            # Get current memory usage
            torch.cuda.set_device(i)
            allocated = torch.cuda.memory_allocated(i) / (1024**3)
            cached = torch.cuda.memory_reserved(i) / (1024**3)
            free = total_mem - cached
            
            gpu_info.append({
                'id': i,
                'name': props.name,
                'total_memory_gb': total_mem,
                'free_memory_gb': free,
                'allocated_gb': allocated,
                'cached_gb': cached,
                'compute_capability': f"{props.major}.{props.minor}"
            })
    return gpu_info

def select_best_gpu():
    """Automatically select the best available GPU"""
    gpu_info = get_gpu_memory_info()
    if not gpu_info:
        return None, "No GPUs available"
    
    # Sort by free memory (descending)
    best_gpu = max(gpu_info, key=lambda x: x['free_memory_gb'])
    return best_gpu['id'], f"Selected GPU {best_gpu['id']} ({best_gpu['name']}) with {best_gpu['free_memory_gb']:.1f}GB free"

def setup_device(gpu_arg):
    """Setup device based on user specification"""
    if gpu_arg.lower() == 'cpu':
        return torch.device('cpu'), "Using CPU (forced by user)"
    
    if not torch.cuda.is_available():
        return torch.device('cpu'), "No GPUs available, falling back to CPU"
    
    if gpu_arg.lower() == 'auto':
        gpu_id, msg = select_best_gpu()
        if gpu_id is not None:
            device = torch.device(f'cuda:{gpu_id}')
            return device, msg
        else:
            return torch.device('cpu'), msg
    
    # Try to use specific GPU ID
    try:
        gpu_id = int(gpu_arg)
        if gpu_id >= torch.cuda.device_count():
            print(f"⚠️ GPU {gpu_id} not available, falling back to auto-selection")
            return setup_device('auto')
        
        device = torch.device(f'cuda:{gpu_id}')
        gpu_info = get_gpu_memory_info()
        gpu_details = next((g for g in gpu_info if g['id'] == gpu_id), None)
        if gpu_details:
            msg = f"Using GPU {gpu_id} ({gpu_details['name']}) with {gpu_details['free_memory_gb']:.1f}GB free"
        else:
            msg = f"Using GPU {gpu_id}"
        return device, msg
    except ValueError:
        print(f"⚠️ Invalid GPU specification '{gpu_arg}', falling back to auto-selection")
        return setup_device('auto')

def optimize_gpu_settings(device, max_memory_fraction=0.8):
    """Optimize GPU settings for memory efficiency"""
    if device.type == 'cuda':
        # Set memory fraction
        if hasattr(torch.cuda, 'set_per_process_memory_fraction'):
            torch.cuda.set_per_process_memory_fraction(max_memory_fraction, device.index)
        
        # Enable memory efficient attention if available
        if hasattr(torch.backends.cuda, 'enable_flash_sdp'):
            torch.backends.cuda.enable_flash_sdp(True)
        
        # Set CUDA memory management
        os.environ.setdefault('PYTORCH_CUDA_ALLOC_CONF', 'max_split_size_mb:512,roundup_power2_divisions:16')

def calculate_optimal_batch_size(device, base_batch_size=16):
    """Calculate optimal batch size based on available GPU memory"""
    if device.type == 'cpu':
        return min(base_batch_size, 8)  # Conservative for CPU
    
    gpu_info = get_gpu_memory_info()
    if not gpu_info:
        return base_batch_size
    
    gpu_details = next((g for g in gpu_info if g['id'] == device.index), None)
    if not gpu_details:
        return base_batch_size
    
    free_memory_gb = gpu_details['free_memory_gb']
    
    # Rough estimation: adjust batch size based on available memory
    if free_memory_gb > 20:
        return base_batch_size
    elif free_memory_gb > 10:
        return base_batch_size
    elif free_memory_gb > 5:
        return max(base_batch_size // 2, 4)
    else:
        return max(base_batch_size // 4, 2)

# ─────────────────────────────────────────────────────────────
# Config
# ─────────────────────────────────────────────────────────────
base_path = "/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq"
base_model_nms = ["TSp_vs_nonProm", "TSp_vs_genNullseqs"]
length_nm = "3k"
subsets = ["tspAll_", "tspliver_", "tsptestis_", "tspbrain_"]

data_path = "/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/data/jul_2025/split"
tokenizer_name = "jaandoui/DNABERT2-AttentionExtracted"

# ─────────────────────────────────────────────────────────────
# Device and Memory Setup
# ─────────────────────────────────────────────────────────────
args = parse_args()
device, device_msg = setup_device(args.gpu)
print(f"🔧 {device_msg}")

# Show all available GPUs for reference
if torch.cuda.is_available():
    print(f"📊 Available GPUs:")
    for gpu in get_gpu_memory_info():
        print(f"   GPU {gpu['id']}: {gpu['name']} ({gpu['total_memory_gb']:.1f}GB total, {gpu['free_memory_gb']:.1f}GB free)")

# Optimize GPU settings
optimize_gpu_settings(device, args.max_memory_fraction)

# Calculate optimal batch size
BATCH_SIZE = calculate_optimal_batch_size(device, args.batch_size)
print(f"🎯 Using batch size: {BATCH_SIZE}")

# ─────────────────────────────────────────────────────────────
# Job Collection (unchanged)
# ─────────────────────────────────────────────────────────────
jobs = []

for base_model_nm in base_model_nms:
    length_dir = os.path.join(base_path, base_model_nm, length_nm)
    if not os.path.exists(length_dir):
        print(f"⚠️ Missing length dir: {length_dir}")
        continue

    all_subdirs = [d for d in os.listdir(length_dir) if os.path.isdir(os.path.join(length_dir, d))]
    matched_subdirs = [d for d in all_subdirs if any(d.startswith(pref) for pref in subsets)]
    if not matched_subdirs:
        print(f"⚠️ No matching subset dirs under {length_dir}")
        continue

    for subset in matched_subdirs:
        search_dir = os.path.join(length_dir, subset)
        lr_dirs = [d for d in os.listdir(search_dir) if d.startswith("lr3e-5")]
        if not lr_dirs:
            print(f"⚠️ No lr3e-5 dirs under {search_dir}")
            continue

        for lr_dir in lr_dirs:
            model_dir = os.path.join(search_dir, lr_dir)
            checkpoints = [d for d in os.listdir(model_dir) if d.startswith("checkpoint-")]
            if not checkpoints:
                print(f"⚠️ No checkpoints in {model_dir}")
                continue

            latest_ckpt = sorted(checkpoints, key=lambda x: int(x.split("-")[1]))[-1]
            model_path = os.path.join(model_dir, latest_ckpt)
            data_dir = os.path.join(data_path, base_model_nm, length_nm, subset)

            res_pdir = f"{base_path}/RESULT/{lr_dir}/{base_model_nm}_{length_nm}_{subset}"
            os.makedirs(res_pdir, exist_ok=True)

            jobs.append({
                "base_model_nm": base_model_nm,
                "subset": subset,
                "lr_dir": lr_dir,
                "model_path": model_path,
                "data_dir": data_dir,
                "res_pdir": res_pdir,
            })

print(f"-> Collected {len(jobs)} jobs")
for j in jobs:
    print(f"** {j['base_model_nm']} | {j['subset']} | {j['lr_dir']} | {os.path.basename(j['model_path'])}")

# ─────────────────────────────────────────────────────────────
# Enhanced Utility Functions with Better Memory Management
# ─────────────────────────────────────────────────────────────
def load_model_and_tokenizer(model_path, tokenizer_name):
    """Load model and tokenizer with memory optimization"""
    print(f"🔄 Loading model from {model_path}")
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_name, trust_remote_code=True)
    
    # Load model with memory optimization
    model = AutoModelForSequenceClassification.from_pretrained(
        model_path, 
        output_attentions=True, 
        trust_remote_code=True,
        torch_dtype=torch.float32, 
        low_cpu_mem_usage=True
    )
    
    model = model.to(device)
    model.eval()
    
    # Print model memory usage
    if device.type == 'cuda':
        allocated = torch.cuda.memory_allocated(device) / (1024**3)
        print(f"📊 Model loaded. GPU memory allocated: {allocated:.2f}GB")
    
    return tokenizer, model

def clear_gpu_cache():
    """Clear GPU cache and run garbage collection"""
    if device.type == 'cuda':
        torch.cuda.empty_cache()
        torch.cuda.ipc_collect()
    gc.collect()

def extract_attention_records(model, tokenizer, data_dir, split_file, batch_size=BATCH_SIZE):
    split = split_file.replace(".csv", "")
    path = os.path.join(data_dir, split_file)
    if not os.path.exists(path):
        print(f"⚠️ Missing file: {path}")
        return []

    df_raw = pd.read_csv(path)
    sequences = df_raw["Sequence"].tolist()
    labels = df_raw["Label"].tolist()

    records = []
    total_batches = (len(sequences) + batch_size - 1) // batch_size
    
    print(f"🔄 Processing {len(sequences)} sequences in {total_batches} batches of size {batch_size}")
    
    for batch_idx, start in enumerate(range(0, len(sequences), batch_size)):
        if batch_idx % 10 == 0:  # Progress update every 10 batches
            print(f"   Batch {batch_idx + 1}/{total_batches}")
            
        batch_seqs = sequences[start:start+batch_size]
        batch_labels = labels[start:start+batch_size]

        try:
            inputs = tokenizer(
                batch_seqs, 
                padding=True, 
                truncation=True, 
                max_length=256, 
                return_tensors="pt"
            ).to(device)
            
            with torch.no_grad():
                outputs = model(**inputs)

            # Move to CPU immediately to save GPU memory
            all_attentions = torch.stack(outputs.attentions, dim=0)  # (layers, batch, heads, seq, seq)
            all_attentions = all_attentions.cpu()  # Move to CPU
            num_layers = all_attentions.shape[0]

            for layer_idx in range(num_layers):
                avg_attn = all_attentions[layer_idx].mean(dim=1)  # (batch, seq, seq)
                for i, seq in enumerate(batch_seqs):
                    token_ids = inputs["input_ids"][i].cpu()
                    tokens = tokenizer.convert_ids_to_tokens(token_ids)
                    attn_score = avg_attn[i].sum(dim=0)
                    attn_score = attn_score / attn_score.max() if attn_score.max() > 0 else attn_score
                    attn_score_list = attn_score.numpy().tolist()

                    records.append({
                        "split": split,
                        "layer": layer_idx,
                        "label": batch_labels[i],
                        "sequence": seq,
                        "tokens": tokens,
                        "attn_scores": attn_score_list,
                    })

            # Cleanup after each batch
            del inputs, outputs, all_attentions
            clear_gpu_cache()
            
        except torch.cuda.OutOfMemoryError as e:
            print(f"❌ GPU OOM error in batch {batch_idx + 1}. Try reducing batch size.")
            print(f"   Current GPU memory: {torch.cuda.memory_allocated(device) / (1024**3):.2f}GB")
            raise e
        except Exception as e:
            print(f"❌ Error in batch {batch_idx + 1}: {e}")
            continue

    return records

def expand_token_records(row):
    try:
        tokens = ast.literal_eval(str(row["tokens"]))
        scores = ast.literal_eval(str(row["attn_scores"]))
    except (ValueError, SyntaxError):
        return pd.DataFrame()
    if len(tokens) != len(scores):
        return pd.DataFrame()

    return pd.DataFrame({
        "split": [row["split"]] * len(tokens),
        "layer": [row["layer"]] * len(tokens),
        "label": [row["label"]] * len(tokens),
        "sequence": [row["sequence"]] * len(tokens),
        "token": tokens,
        "attn_score": scores,
    })

def compute_token_bounds(row):
    sequence, token = row["sequence"], row["token"]
    start = sequence.find(token)
    end = start + len(token) if start != -1 else None
    return pd.Series({"start_pos": start + 1, "end_pos": end + 1})

def apply_attention_mask(group):
    scores = group["attn_score"]
    token_len = group["token"].astype(str).apply(len)
    mean_attn, min_attn = scores.mean(), scores.min()
    mask = (scores > mean_attn) & (scores > 2 * min_attn) & (token_len >= 4)
    group = group.copy()
    group["mask"] = mask.astype(int)
    return group

def process_job(job, tokenizer_name):
    raw_path   = os.path.join(job["res_pdir"], "1_Attention_raw.parquet")
    motif_path = os.path.join(job["res_pdir"], "2_Attention_motif.parquet")

    # Skip if already processed
    if os.path.exists(raw_path) and os.path.exists(motif_path):
        print(f"⏭️ Skipping {job['base_model_nm']} | {job['subset']} | {job['lr_dir']} (files already exist)")
        return None, None

    print(f"\n▶️ Processing: {job['base_model_nm']} | {job['subset']} | {job['lr_dir']}")
    
    # Show initial GPU memory
    if device.type == 'cuda':
        initial_memory = torch.cuda.memory_allocated(device) / (1024**3)
        print(f"📊 Initial GPU memory: {initial_memory:.2f}GB")

    tokenizer, model = load_model_and_tokenizer(job["model_path"], tokenizer_name)

    all_layer_records = []
    for split_file in ["train.csv", "dev.csv", "test.csv"]:
        print(f"🔄 Processing {split_file}")
        recs = extract_attention_records(model, tokenizer, job["data_dir"], split_file, batch_size=BATCH_SIZE)
        all_layer_records.extend(recs)
        print(f"✅ {split_file}: {len(recs)} records")

    print(f"🔄 Creating DataFrames...")
    raw_df = pd.DataFrame(all_layer_records)

    print(f"🔄 Expanding token records...")
    motif_raw_df = pd.concat([expand_token_records(row) for _, row in raw_df.iterrows()], ignore_index=True)
    motif_raw_df = motif_raw_df[~motif_raw_df["token"].astype(str).str.startswith("[")]
    
    print(f"🔄 Computing token bounds...")
    motif_raw_df[["start_pos", "end_pos"]] = motif_raw_df.apply(compute_token_bounds, axis=1)
    motif_raw_df["token_length"] = motif_raw_df["end_pos"].astype(int) - motif_raw_df["start_pos"].astype(int)

    print(f"🔄 Applying attention masks...")
    motif_raw_df2 = motif_raw_df.groupby("sequence", group_keys=False).apply(apply_attention_mask)

    print(f"💾 Saving results...")
    raw_df.to_parquet(raw_path, index=False)
    motif_raw_df2.to_parquet(motif_path, index=False)

    print(f"✅ Saved: {raw_path} ({len(raw_df)} records)")
    print(f"✅ Saved: {motif_path} ({len(motif_raw_df2)} records)")

    # Final cleanup
    del model, tokenizer, all_layer_records, raw_df, motif_raw_df, motif_raw_df2
    clear_gpu_cache()
    
    if device.type == 'cuda':
        final_memory = torch.cuda.memory_allocated(device) / (1024**3)
        print(f"📊 Final GPU memory: {final_memory:.2f}GB")

    return None, None  # Don't return large DataFrames to save memory


# ─────────────────────────────────────────────────────────────
# Main Execution
# ─────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print(f"🚀 Starting DNABERT2 attention extraction on {device}")
    print(f"📋 Collected {len(jobs)} jobs with batch size {BATCH_SIZE}")
    
    successful_jobs = 0
    failed_jobs = 0

    # Run only jobs[7] if it exists
    if len(jobs) > 5:
        job = jobs[5]
        try:
            print(f"\n{'='*60}")
            print(f"🔄 Running single job 7/{len(jobs)}: {job['base_model_nm']} | {job['subset']} | {job['lr_dir']}")
            print(f"{'='*60}")

            process_job(job, tokenizer_name)
            successful_jobs += 1

            if device.type == 'cuda':
                memory_used = torch.cuda.memory_allocated(device) / (1024**3)
                print(f"📊 Current GPU memory: {memory_used:.2f}GB")

        except Exception as e:
            print(f"❌ Error in job 7: {e}")
            failed_jobs += 1
            clear_gpu_cache()
    else:
        print(f"❌ jobs[7] does not exist (only {len(jobs)} jobs collected).")

    print(f"\n{'='*60}")
    print(f"🏁 Processing complete!")
    print(f"✅ Successful jobs: {successful_jobs}")
    print(f"❌ Failed jobs: {failed_jobs}")
    print(f"📊 Total jobs: {len(jobs)}")
    print(f"{'='*60}")
