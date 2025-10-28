#!/bin/bash
set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CONDA_ENV="dnabert_aug_2025_jupyter"

PYTHON_SCRIPT_1="TSProm/src/3_attention/1_raw_attention_extract-gpu.py"
PYTHON_SCRIPT_2A="TSProm/src/3_attention/2A_save_meme.py"
PYTHON_SCRIPT_3="TSProm/src/3_attention/3_SHAP.py"

GPU_ID="auto"
BATCH_SIZE=64

JOBS_FILE="TSProm/src/files/jobs.json"
BASE_PATH="TSProm/data"
DATA_PATH="TSProm/data/split_into_train_test_dev"
LENGTH_NM="3k"
BASE_MODEL_NMS="TSp_vs_nonProm,TSp_vs_genNullseqs"
SUBSETS="tspAll_,tspliver_,tsptestis_,tspbrain_,tspmuscle_"

# need to download from https://jaspar.elixir.no/downloads/
JASPAR_FILE="jaspar/JASPAR2024_CORE_non-redundant_pfms_meme.txt"



echo "Starting DNABERT2 Motif Analysis Pipeline..."


echo "=== STEP 1: Raw Attention Extraction ==="
conda activate "$CONDA_ENV"
python "$PYTHON_SCRIPT_1" --gpu "$GPU_ID" --batch-size "$BATCH_SIZE"


echo "=== STEP 2A: Create Jobs File ==="
conda activate "$CONDA_ENV"
python "$PYTHON_SCRIPT_2A" \
  --jobs_file "$JOBS_FILE" \
  --refresh_jobs \
  --base_path "$BASE_PATH" \
  --data_path "$DATA_PATH" \
  --length_nm "$LENGTH_NM" \
  --base_model_nms "$BASE_MODEL_NMS" \
  --subsets "$SUBSETS" \
  --job_index 0

if [ ! -f "$JOBS_FILE" ]; then
    echo "ERROR: $JOBS_FILE not created."; exit 1
fi

N_JOBS=$(jq 'length' "$JOBS_FILE")
echo "Found $N_JOBS jobs."

for (( i=0; i<$N_JOBS; i++ )); do
    echo "=== Processing Job $i / $((N_JOBS - 1)) ==="
    JOB_RES_DIR=$(jq -r ".[$i].res_pdir" "$JOBS_FILE")

    echo "--- [Job $i] Step 2A (Motif Processing) ---"
    conda activate "$CONDA_ENV_DNABERT"
    python "$PYTHON_SCRIPT_2A" --jobs_file "$JOBS_FILE" --job_index "$i"
    conda deactivate

    echo "--- [Job $i] Step 2B (MEME) ---"
    motif_dir="${JOB_RES_DIR}/motifs"
    fasta_file="${motif_dir}/weblogo_label_1.fa"

    if [ -f "$fasta_file" ]; then
        conda activate "$CONDA_ENV"
        cd "$motif_dir"
        fasta-get-markov weblogo_label_1.fa bg_model_label1.txt
        meme weblogo_label_1.fa -dna -oc meme_out_label_1 -mod zoops -nmotifs 50 -minw 6 -maxw 20 -revcomp -bfile bg_model_label1.txt
        tomtom -oc tomtom_out_label_1 -evalue -thresh 1.0 meme_out_label_1/meme.txt "$JASPAR_FILE"
        cd "$SCRIPT_DIR"
    else
        echo "Skipping MEME: FASTA not found ($fasta_file)"
    fi

    echo "--- [Job $i] Step 3 (SHAP) ---"
    shap_input_csv="${JOB_RES_DIR}/4sig_motifs_prePWM.csv"
    if [ -f "$shap_input_csv" ]; then
        conda activate "$CONDA_ENV_DNABERT"
        TEMP_SHAP_SCRIPT="${SCRIPT_DIR}/temp_shap_job_${i}.py"
        sed "s|jobs_path = \".*\"|jobs_path = \"$JOBS_FILE\"|" "$PYTHON_SCRIPT_3" > "$TEMP_SHAP_SCRIPT"
        sed -i "s|load_job_from_json(jobs_path, index=.*)|load_job_from_json(jobs_path, index=$i)|" "$TEMP_SHAP_SCRIPT"
        python "$TEMP_SHAP_SCRIPT"
        rm "$TEMP_SHAP_SCRIPT"
        conda deactivate
    else
        echo "Skipping SHAP: CSV not found ($shap_input_csv)"
    fi
done

echo "Pipeline Complete."
