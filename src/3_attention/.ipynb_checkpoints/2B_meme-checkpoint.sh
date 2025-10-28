# conda activate scrna-seq
# bash /data/private/psurana/TSpDNA2/src/3_attention/2B_meme.sh

res_pdir="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq/RESULT/lr3e-5_ep10/"
dirs=(
  "TSp_vs_nonProm_3k_tspAll_nonPromHu"
  "TSp_vs_nonProm_3k_tspliver_nonPromHu"
  "TSp_vs_genNullseqs_3k_tspAll_genNullseqs"
  "TSp_vs_genNullseqs_3k_tspliver_genNullseqs"
  "TSp_vs_genNullseqs_3k_tsptestis_genNullseqs"
)
dirs=(
  "TSp_vs_nonProm_3k_tspmuscle_nonPromHu"
  "TSp_vs_genNullseqs_3k_tspmuscle_genNullseqs"
)
for dir in "${dirs[@]}"; do
  cd "$res_pdir/$dir/motifs"

  # background model
  fasta-get-markov weblogo_label_1.fa bg_model_label1.txt

  # meme for novel motif discovery
  meme weblogo_label_1.fa \
    -dna \
    -oc meme_out_label_1 \
    -mod zoops \
    -nmotifs 50 \
    -minw 6 \
    -maxw 20 \
    -revcomp \
    -bfile bg_model_label1.txt

  # tomtom
  tomtom -oc tomtom_out_label_1 -evalue -thresh 1.0 \
    meme_out_label_1/meme.txt \
    /data/projects/dna/pallavi/jaspar/JASPAR2024_CORE_non-redundant_pfms_meme.txt
done


# run tomtom --> compares with JASPAR2024
# cd /data/projects/dna/pallavi/jaspar
# wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_meme.zip
# unzip JASPAR2024_CORE_non-redundant_pfms_meme.zip -d JASPAR2024_MEME

# wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_meme.txt
