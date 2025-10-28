

# --model_name_or_path zhihan1996/DNABERT-2-117M \
# --model_name_or_path jaandoui/DNABERT2-AttentionExtracted \


# transformers==4.29.2


##-------------------------------------------------------------------------------------------------------------------
# try diff lengths for specific folders on both models
##-------------------------------------------------------------------------------------------------------------------
set -euo pipefail 

cd /data/projects/dna/pallavi/DNABERT_runs

export CUDA_VISIBLE_DEVICES=5,6,7
export NEPTUNE_PROJECT="ladfgbdf/sfsdf"
# export NEPTUNE_PROJECT="lavi.ana/genNUllseqs"

#--------------------------------------------------
# Aug 26 2025
#--------------------------------------------------


# conda activate dnabert_aug_2025
export TOKENIZERS_PARALLELISM=false
export NEPTUNE_API_TOKEN="eyQ=="

# Base data path
base_path="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/data/jul_2025/TSp_vs_nonProm"

# List of datasets (all tissue-type directories)
datasets=("tsp_spleen_null" "tsp_spleen_wide" "tsp_spleen_tenh" "tsp_spleen_low" "tsp_muscle_null" "tsp_muscle_wide" "tsp_muscle_tenh" "tsp_muscle_low")

# Ordered list of lengths (to process same lengths together across datasets)
lengths=("2000" "3000" "4000")

# Hyperparameters
learning_rates=(3e-6 5e-6 7e-6 1e-5 3e-5 1e-4 3e-4)

declare -p learning_rates

default_batch_size=8
testis_batch_size=4
gradient_accumulation_steps=1
num_train_epochs=10
warmup_steps=10
weight_decay=1e-2

: "${weight_decay:=1e-2}"

##------------------------------------------------------------##

# Length-to-MaxLength Mapping (Modify if needed)
declare -A length_map
length_map=( ["2000"]=500 ["3000"]=750 ["4000"]=1000 ["550"]=138 ["750"]=188 ["350"]=88 ["1200"]=300 )


# Process each length across all datasets first
for length in "${lengths[@]}"; do
    length="${length%,}"
    echo "Processing all datasets for length: $length"
    
    # Loop through each dataset (tissue-type directory)
    for data in "${datasets[@]}"; do
        dataset_path="${base_path}/${length}/${data}"

        # Check if dataset-length directory exists
        if [ ! -d "$dataset_path" ]; then
            echo "Skipping: $dataset_path (does not exist)"
            continue
        fi
        
        # Choose batch size based on tissue type
        if [[ "$data" == *"testis"* ]]; then
            batch_size=$testis_batch_size
        else
            batch_size=$default_batch_size
        fi

        # Set MAX_LENGTH dynamically
        MAX_LENGTH=${length_map[$length]}
        
        # Loop through learning rates
        for lr in "${learning_rates[@]}"; do
            run_name="DNA2_${data}_${length}_lr${lr}_ep${num_train_epochs}"
            output_dir="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/unique_TSS/human_mouse_mmseq2/${length}/${data}/lr${lr}_ep${num_train_epochs}"
            
            # Check if eval_results.json already exists
            results_file="${output_dir}/results/${run_name}/eval_results.json"
            if [ -f "$results_file" ]; then
                echo "✅ Skipping ${data}/${length}/lr${lr} — results already exist: $results_file"
                continue
            fi

            # Create output directory
            mkdir -p "$output_dir"

            echo "🚀 Running finetuning for ${data}/${length} at LR ${lr}"
            python DNABERT2_AttentionExtracted.py \
                --model_name_or_path jaandoui/DNABERT2-AttentionExtracted \
                --data_path "$dataset_path" \
                --kmer -1 \
                --run_name "$run_name" \
                --model_max_length "$MAX_LENGTH" \
                --per_device_train_batch_size "$batch_size" \
                --per_device_eval_batch_size "$batch_size" \
                --gradient_accumulation_steps "$gradient_accumulation_steps" \
                --learning_rate "$lr" \
                --num_train_epochs "$num_train_epochs" \
                --fp16 \
                --save_steps 100 \
                --output_dir "$output_dir" \
                --evaluation_strategy "steps" \
                --eval_steps 100 \
                --warmup_steps "$warmup_steps" \
                --logging_steps 100 \
                --overwrite_output_dir \
                --log_level "info" \
                --weight_decay "$weight_decay"
        done
    done
done


# conda activate dnabert_aug_2025
# cd /data/private/psurana/TSProm/src/1_finetune/DNABERT2/; bash 2_finetune_dna2_human_mouse_MMSEQ2.sh


