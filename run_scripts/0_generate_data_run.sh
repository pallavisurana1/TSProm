#!/bin/bash

# =======================================================
# Bash config for make_data.R
# =======================================================

# ---- CONFIG VARIABLES ----
cd TSProm/src/0_generate_data

export benchmark="make_data.R"
export modelA="make_modelA_negSet.sh"
export modelB="make_modelB_nullSeqs.R"

export SAVE_DIR="/TSProm/data"

export HUMAN_CSV="/vast/projects/rdavuluri-group/psurana/projects/Davuluri_lab/TransTEx_GTEX_Feb_2025/TransTEx_Normal_GTEx_gene_ENSEMBL113.csv"
export MOUSE_RDS="/vast/projects/rdavuluri-group/psurana/projects/Davuluri_lab/Extramapper_2024/TransTEx/MouseTranstex/mouse_transtex_biomart.rds"
export NONPROM_RDS="/vast/projects/rdavuluri-group/psurana/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/rds/hg38_nonPromoter_randomSeqs.rds"

export SPECIES="human,mouse"
export TISSUES="brain,testis,liver"
export TRAIN_P=0.8
export DEV_P=0.2
export TSS_MIN_SEP=50


# ---- RUN PIPELINE for making benchmark data of tsp vs rest -----
Rscript "$benchmark"

# ---- RUN PIPELINE for making neg set for model A -----
bash "$modelA" \
  --workdir /data/projects/promoters/raw1 \
  --outdir /data/projects/promoters/output \
  --species human \
  --windows w1k1k,w3k1k

# ---- RUN PIPELINE for making neg set for model B -----
Rscript "$modelB" \
  --input_rds /path/to/human_seq.rds \
  --outdir /path/to/out/run1 \
  --tissues spleen,muscle \
  --windows 1k,2k,3k