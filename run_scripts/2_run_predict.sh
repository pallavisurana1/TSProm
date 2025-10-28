set -euo pipefail

conda activate dnabert_aug_2025
export CUDA_VISIBLE_DEVICES=5,6,7
export TOKENIZERS_PARALLELISM=false

base_model_nm="TSp_vs_genNullseqs"
subset="muscle_genNullseqs"
lr_dir="lr3e-5_ep10"
model_path="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq/TSp_vs_genNullseqs/3k/muscle_genNullseqs/lr3e-5_ep10/checkpoint-300"
data_dir="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/data/jul_2025/split/TSp_vs_genNullseqs/3k/muscle_genNullseqs"
res_pdir="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq/RESULT/lr3e-5_ep10/TSp_vs_genNullseqs_3k_tspmuscle_genNullseqs"

python TSProm/src/2_predict/1_predict.py \
  --model_name_or_path $model_path \
  --data_path $data_dir \
  --output_dir $res_pdir/preds \
  --per_device_eval_batch_size 1 \
  --model_max_length 750 \
  --kmer -1


# to run multiple predictions
python TSProm/src/2_predict/1_runAll_predict.py