#!/usr/bin/env python3
import json
import subprocess

with open("/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq/jobs.json", "r") as f:
    jobs = json.load(f)

script_path = "/data/private/psurana/DNA2-TSp/code/predict/1_predict.py"

for job in jobs:
    model_path = job["model_path"]
    data_dir = job["data_dir"]
    output_dir = f'{job["res_pdir"]}/preds'

    cmd = [
        "python", script_path,
        "--model_name_or_path", model_path,
        "--data_path", data_dir,
        "--output_dir", output_dir,
        "--per_device_eval_batch_size", "1",
        "--model_max_length", "750",
        "--kmer", "-1"
    ]

    print(f"🔁 Running prediction for {job['subset']}...")
    subprocess.run(cmd)
    print(f"✅ Done: {job['subset']}\n")


# conda activate dnabert_aug_2025; export CUDA_VISIBLE_DEVICES=1,2,3,4,5,6,7; export TOKENIZERS_PARALLELISM=false; python /data/private/psurana/DNA2-TSp/code/predict/1_runAll_predict.py