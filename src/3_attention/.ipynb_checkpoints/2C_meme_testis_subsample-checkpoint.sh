# conda activate scrna-seq
# bash /data/private/psurana/TSpDNA2/src/3_attention/2B_meme.sh

#!/bin/bash

res_pdir="/data/projects/dna/pallavi/DNABERT_runs/DATA_RUN/dnabert2_FineTune_Zhihan_attention_extracted/july_2025_mmseq/RESULT/lr3e-5_ep10/"
dir="TSp_vs_nonProm_3k_tsptestis_nonPromHu"
full_fasta_file="$res_pdir/$dir/motifs/weblogo_label_1.fa"
jaspar_db="/data/projects/dna/pallavi/jaspar/JASPAR2024_CORE_non-redundant_pfms_meme.txt"

# --- ITERATIVE SUBSAMPLING LOOP ---
# Let's run it 10 times
for i in {1..10}
do
  echo "--- Starting Run $i ---"
  
  # Define a unique output directory for this run
  run_dir="$res_pdir/$dir/motifs/subsample_run_$i"
  mkdir -p "$run_dir"
  cd "$run_dir"

  # 1. Create a random subsample of 400 sequences
  # The -s flag sets a random seed, making each run different
  subsample_fa="subsample_$i.fa"
  /home/campus.stonybrook.edu/psurana/seqtk/seqtk sample -s$i "$full_fasta_file" 1000 > "$subsample_fa"

  # 2. Create the background model from the subsample
  fasta-get-markov "$subsample_fa" bg_model.txt

  # 3. Run MEME on the subsample
  meme "$subsample_fa" \
    -dna \
    -oc meme_out \
    -mod zoops \
    -nmotifs 10 \
    -minw 6 \
    -maxw 20 \
    -revcomp \
    -bfile bg_model.txt

  # 4. Run Tomtom on the results of the subsample
  tomtom -oc tomtom_out -evalue -thresh 1.0 \
    meme_out/meme.txt \
    "$jaspar_db"
  
  echo "--- Finished Run $i ---"
done

echo "All subsampling runs complete. Now aggregate and analyze the results from the 'subsample_run_*' directories."

# run tomtom --> compares with JASPAR2024
# cd /data/projects/dna/pallavi/jaspar
# wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_meme.zip
# unzip JASPAR2024_CORE_non-redundant_pfms_meme.zip -d JASPAR2024_MEME

# wget https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_meme.txt
