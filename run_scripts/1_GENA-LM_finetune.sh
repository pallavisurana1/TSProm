#!/usr/bin/env bash
set -euo pipefail

#############################################
#  Run GenaLM fine-tuning job (multi-GPU)
#  Usage:
#    bash run_finetune_genalm.sh
#  or for background run with log:
#    nohup bash run_finetune_genalm.sh > finetune.log 2>&1 &
#############################################

# ==== 1. CONFIGURATION ====
CONDA_ENV="dnabert2.0"
PROJECT_DIR="TSProm/src/1_finetune"
SCRIPT="FineTune_GENALM.py"

# Neptune setup
export NEPTUNE_PROJECT="fbghfgb.sdvsdv/sdvsdv-human-mouse-Finetune"
export NEPTUNE_API_TOKEN="F="

# GPU & Python runtime
export CUDA_VISIBLE_DEVICES="5,6,7"
export TOKENIZERS_PARALLELISM="false"
export PYTORCH_CUDA_ALLOC_CONF="expandable_segments:True"

# ==== 2. ACTIVATE ENVIRONMENT ====
echo ">>> Activating conda env: $CONDA_ENV"
source ~/.bashrc
conda activate "$CONDA_ENV"

# ==== 3. RUN TRAINING ====
echo ">>> Starting fine-tuning at $(date)"
cd "$PROJECT_DIR"

python "$SCRIPT"

echo ">>> Training finished at $(date)"
