# conda activate R411

# Activate the environment
# conda activate R411

# Load necessary libraries
pacman::p_load(data.table, stringr, dplyr, readxl, feather, biomaRt, tidyverse, openxlsx, ggplot2, VennDiagram, janitor, gridExtra)
library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmusculus.UCSC.mm39)  # Mouse genome
library(biomaRt)

# Try the US Ensembl mirror
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                   host = "https://useast.ensembl.org")

# Load the human dataset
human_ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
mouse_ensembl <- useDataset("mmusculus_gene_ensembl", mart = ensembl)

# Define attributes to retrieve
attributes <- c("chromosome_name", "ensembl_transcript_id", 
                "transcription_start_site", "strand", "external_gene_name")

# Download human data
human_data <- getBM(attributes = attributes, mart = human_ensembl)

# Download mouse data
mouse_data <- getBM(attributes = attributes, mart = mouse_ensembl)

# Filter valid chromosomes (autosomes + X, Y)
valid_chromosomes <- c(as.character(1:22), "X", "Y")

human_data <- human_data %>% filter(chromosome_name %in% valid_chromosomes)
mouse_data <- mouse_data %>% filter(chromosome_name %in% as.character(1:19) | chromosome_name %in% c("X", "Y"))

# Format and calculate coordinates
process_data <- function(df) {
  df %>%
    mutate(
      chromosome_name = str_c("chr", chromosome_name),
      strand = ifelse(strand == 1, "+", "-"),
      start_modified = ifelse(strand == "+", transcription_start_site - 5000, transcription_start_site - 999),
      end_modified = ifelse(strand == "+", transcription_start_site + 999, transcription_start_site + 5000)
    ) %>%
    mutate(
      detail = paste(chromosome_name, start_modified, end_modified,
                     ensembl_transcript_id, transcription_start_site, strand,
                     sep = "_")
    ) %>%
    select(chromosome_name, start_modified, end_modified,
           ensembl_transcript_id, strand, transcription_start_site, detail)
}

human_processed <- process_data(human_data)
mouse_processed <- process_data(mouse_data)

# Save as BED files
human_bed_path <- "/home/psurana/projects/Davuluri_lab/DNABERT_TSpProm/Nov_2024/dnabert2_finetune_ensembl113/5k_TSS_1k/human/raw.bed"
mouse_bed_path <- "/home/psurana/projects/Davuluri_lab/DNABERT_TSpProm/Nov_2024/dnabert2_finetune_ensembl113/5k_TSS_1k/mouse/raw.bed"

write_tsv(human_processed, human_bed_path, col_names = FALSE)
write_tsv(mouse_processed, mouse_bed_path, col_names = FALSE)


##--------------------------------------------------------------------------------------------------------
## Convert BED to FASTA for both human and mouse
##--------------------------------------------------------------------------------------------------------

convert_bed_to_fasta <- function(bed_path, genome, output_fasta) {
  # Read the BED file
  bed_data <- read.table(bed_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(bed_data) <- c("chromosome_name", "start_modified", "end_modified",
                          "transcript", "strand", "transcription_start_site", "detail")
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = bed_data$chromosome_name,
    ranges = IRanges(start = bed_data$start_modified, end = bed_data$end_modified),
    strand = bed_data$strand
  )
  
  # Fetch sequences
  sequences <- getSeq(genome, gr)
  
  # Create headers for FASTA
  headers <- apply(bed_data[c("chromosome_name", "start_modified", "end_modified",
                              "transcript", "transcription_start_site", "strand")], 1, function(row) {
                                paste(row, collapse = "|")
                              })
  
  # Remove spaces
  headers <- gsub(" ", "", headers)
  names(sequences) <- headers
  
  # Save FASTA
  writeXStringSet(sequences, filepath = output_fasta)
  message("Saved FASTA: ", output_fasta)
}

# Convert to FASTA
convert_bed_to_fasta(human_bed_path, BSgenome.Hsapiens.UCSC.hg38, "6k.fasta")
convert_bed_to_fasta(mouse_bed_path, BSgenome.Mmusculus.UCSC.mm39, "6k.fasta")

