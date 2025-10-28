
#-conda activate R411

pacman::p_load(data.table, stringr, dplyr,readxl, feather, biomaRt, tidyverse, openxlsx, ggplot2, VennDiagram, janitor, gridExtra)
library(GenomicRanges)
library(Biostrings)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(caret)

##--------------------------------------------------------------------------------------------------------
## make bed file 
##--------------------------------------------------------------------------------------------------------

# 1k_TSS_1k  2k_TSS_1k  3k_TSS_1k  500_TSS_50

##-paths
ft = "/home/psurana"
# ft = "Z:"
res_path =paste0(ft, "/projects/Davuluri_lab/2_Tissue-specific-July-2022/subsample/manuscript/")
save_pth_dna = paste0(ft, "/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/")

##--------------------------------------------------------------------------------------------------------

# red human and mouse transtex
df_human=fread(file.path(ft, "/projects/Davuluri_lab/TransTEx_GTEX_Feb_2025/TransTEx_Normal_GTEx_gene_ENSEMBL113.csv"))
df_mouse=fread(file.path(ft, "/projects/Davuluri_lab/CancerTransTEx/RES/solid_tumors/cancer_transtex_gene_ENSEMBL113.csv"))
dim(df_human); dim(df_mouse)

##--------------------------------------------------------------------------------------------------------

# remove chr mt and redundant cols
df_human <- df_human %>% filter(chromosome_name %in% as.character(1:22) | chromosome_name %in% c("X", "Y"))
df_mouse <- df_mouse %>% filter(chromosome_name %in% as.character(1:19) | chromosome_name %in% c("X", "Y"))

df = df_human

# calculate TSS differences
df1 <- df %>%
  arrange(`transcription_start_site`) %>% # Arrange in ascending order
  mutate(TSS_diff = c(0, diff(`transcription_start_site`))) # Calculate differences

# Extract TSS difference and clean data
df2 <- df1 %>%
  mutate(
    chromosome_name = paste0("chr", chromosome_name),
    start = as.integer(transcript_start),
    end = as.integer(transcript_end),
    transcription_start_site = as.integer(transcription_start_site)
  ) %>%
  arrange(TSS_diff) %>%
  distinct(transcription_start_site, .keep_all = TRUE)

df2 = df2 %>% mutate(strand = ifelse(strand == 1, "+", "-"))

# make start and end -10k TSS +2k
data_bed <- df2 %>%
  mutate(start_modified = if_else(strand == "+",       
                                  transcription_start_site - 10000, transcription_start_site - 1999),
         end_modified = if_else(strand == "+", 
                                transcription_start_site + 1999, transcription_start_site + 10000))
  
#  save tsp tissue wise
data_bed_tsp = data_bed %>% filter(category_final == "TspTrans"); dim(data_bed_tsp) 
data_bed_tsp1 = data_bed_tsp %>% dplyr::select(chromosome_name, start_modified, end_modified, 
                                        transcript,transcription_start_site, strand,
                                        tissue, ensembl_gene_id, TSS_diff) 

# colnames for this bed file
data_bed_colnames <- c("chromosome_name", "start_modified", "end_modified", 
                            "transcript", "transcription_start_site", "strand", 
                            "tissue", "ensembl_gene_id", "TSS_diff")


# save bed
bed_nm = "10k_TSS_2k_Correct" #5000 1k_1k
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Lung")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Lung.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Pancreas")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Pancreas.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Spleen")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Spleen.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Nerve")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Nerve.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Thyroid")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Thyroid.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Muscle")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Muscle.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Pituitary")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Pituitary.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Brain")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Brain.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Liver")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Liver.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
write.table(na.omit(subset(data_bed_tsp1, data_bed_tsp1$tissue == "Testis")), file = paste0(save_pth_dna, "bed/", bed_nm, "/", "TspTrans_Testis.bed"), quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)


# save expression groups 
categories <- c("TspTrans", "TenhTrans", "WideTrans", "NullTrans", "LowTrans")
for (category in categories) {
  data_bed_filtered <- data_bed %>% filter(category_final == category)
  data_bed_selected <- data_bed_filtered %>% dplyr::select(all_of(data_bed_colnames))
  
  write.table(data_bed_selected, 
              file = paste0(save_pth_dna, "bed/", bed_nm, "/", category, "_all.bed"), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

##--------------------------------------------------------------------------------------------------------
## bed to fasta
##--------------------------------------------------------------------------------------------------------

bed_dir <- paste0(save_pth_dna, "bed/", bed_nm)
output_fasta_dir <- paste0(save_pth_dna, "fasta/", bed_nm, "/")
dir.create(output_fasta_dir, recursive = TRUE, showWarnings = FALSE)

# Define the reference genome (BSgenome is faster if available)
genome <- BSgenome.Hsapiens.UCSC.hg38

# Function to read a BED file, convert to GRanges, fetch FASTA, and save
convert_bed_to_fasta <- function(bed_file, genome, output_dir) {
  # Read the BED file
  bed_data <- read.table(bed_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  data_bed_colnames <- c("chromosome_name", "start_modified", "end_modified", 
                                               "transcript", "transcription_start_site", "strand", 
                                               "tissue", "ensembl_gene_id", "TSS_diff")
  colnames(bed_data) <- data_bed_colnames
  
  # Convert to GRanges
  gr <- GRanges(
    seqnames = bed_data$chromosome_name,
    ranges = IRanges(start = bed_data$start_modified, end = bed_data$end_modified),
    strand = bed_data$strand
  )
  
  # Fetch sequences
  sequences <- getSeq(genome, gr)
  
  # Select the columns for the header
  header_columns <- c(
    "chromosome_name", "start_modified", "end_modified",
    "transcript", "transcription_start_site", "strand",
    "ensembl_gene_id"
  )
  
  # Create headers by collapsing the selected columns with '|'
  headers <- apply(bed_data[header_columns], 1, function(row) {
    paste(row, collapse = "|")
  })
  
  # Remove any unintended spaces in the headers
  headers <- gsub(" ", "", headers);
  names(sequences) <- headers

  # Define output file path
  fasta_file <- paste0(output_fasta_dir, basename(tools::file_path_sans_ext(bed_file)), ".fasta")
  
  # Save sequences to FASTA
  writeXStringSet(sequences, filepath = fasta_file)
  message("Saved FASTA: ", fasta_file)
}

# Process all BED files
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE); bed_files
lapply(bed_files, function(bed_file) convert_bed_to_fasta(bed_file, genome, output_fasta_dir))



##--------------------------------------------------------------------------------------------------------
## DNABERt 2.0 data - need sequence and label
##--------------------------------------------------------------------------------------------------------

bed_nm = "10k_TSS_2k_Correct" 

# Paths
output_fasta_dir <- paste0(save_pth_dna, "fasta/", bed_nm, "/")

# Function to read FASTA and create dataframe
read_fasta_as_df <- function(fasta_path) {
  fasta_data <- readDNAStringSet(fasta_path)
  data.frame(
    Name = names(fasta_data),
    Sequence = as.character(fasta_data),
    stringsAsFactors = FALSE
  )
}

# List all FASTA files in the directory
fasta_files <- list.files(output_fasta_dir, pattern = "\\.fasta$", full.names = TRUE); fasta_files

# Read all FASTA files into a list of dataframes
fasta_dfs <- lapply(fasta_files, function(fasta_path) {
  read_fasta_as_df(fasta_path)
})

# Extract TSS from the `Name` column and create a new column `TSS`
fasta_dfs <- lapply(fasta_dfs, function(df) {
  df$TSS <- sapply(strsplit(df$Name, "\\|"), function(x) as.numeric(x[5]))
  return(df)
}); names(fasta_dfs) <- sub("\\.fasta$", "", basename(fasta_files))


# Define groups for Label = 0
label_0 <- c("LowTrans_all", "NullTrans_all", "TenhTrans_all", "WideTrans_all")

# Add labels to dataframes in fasta_dfs
fasta_dfs <- lapply(names(fasta_dfs), function(name) {
  fasta_dfs[[name]]$Label <- ifelse(name %in% label_0, 0, 1)
  fasta_dfs[[name]]
}); names(fasta_dfs) <- sub("\\.fasta$", "", basename(fasta_files))


##-------------------
# combine label 0 and 1
##-------------------

# Define groups to combine
label_0 <- c("LowTrans_all", "NullTrans_all", "TenhTrans_all", "WideTrans_all")
tsp_groups <- setdiff(names(fasta_dfs), label_0) 

combined_fasta_dfs <- list()

for (base in label_0) {
  for (tsp in tsp_groups) {
   
    combined_name <- paste0(base, "_", tsp)
    combined_fasta_dfs[[combined_name]] <- rbind(
      fasta_dfs[[base]],
      fasta_dfs[[tsp]]
    )}
  }

##-------------------
# Function to calculate TSS differences
##-------------------

calculate_tss_difference <- function(df, threshold = 50) {
  df <- df[order(df$TSS), ]  # Sort by TSS
  df$TSS_diff <- c(NA, diff(df$TSS))  # Calculate differences
  df <- df[is.na(df$TSS_diff) | df$TSS_diff >= threshold, ]  # Filter by threshold
  return(df)
}

# Function to balance rows with Label 0 and Label 1
balance_labels <- function(df) {
  label_0 <- df[df$Label == 0, ]
  label_1 <- df[df$Label == 1, ]
  
  n0 <- nrow(label_0)
  n1 <- nrow(label_1)
  
  if (n0 > n1) {
    label_0 <- label_0[sample(n0, n1), ]  # Downsample label 0
  } else if (n1 > n0) {
    label_1 <- label_1[sample(n1, n0), ]  # Downsample label 1
  }
  
  return(rbind(label_0, label_1))
}

# use func
# Process combined_fasta_dfs to create two lists
list_unique_TSS <- lapply(combined_fasta_dfs, function(df) {
  df <- df[!duplicated(df$TSS), ]  # Keep unique TSS
  df <- df[!grepl("N", df$Sequence), ]  # Remove rows with 'N' in Sequence
  balance_labels(df)  # Balance rows by label
})

list_unique_TSS_TSSdiff50 <- lapply(combined_fasta_dfs, function(df) {
  df <- df[!duplicated(df$TSS), ]  # Keep unique TSS
  df <- df[!grepl("N", df$Sequence), ]  # Remove rows with 'N' in Sequence
  df <- calculate_tss_difference(df, threshold = 50)  # Apply TSS difference
  balance_labels(df)  # Balance rows by label
})

saveRDS(list_unique_TSS_TSSdiff50, paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS_diff_50.rds"))
saveRDS(list_unique_TSS, paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS.rds"))


x=sapply(list_unique_TSS_TSSdiff50, nrow) %>% as.data.frame %>% rownames_to_column("experiment") 
colnames(x)[2]="unique_TSS_diff_50"
y=sapply(list_unique_TSS, nrow) %>% as.data.frame %>% rownames_to_column("experiment")
colnames(y)[2]="unique_TSS"
z=merge(x,y,by="experiment")

##-------------------------------
# train:Dev:test (70:15:15) split
##-------------------------------

# Function to perform stratified sampling and save files with/without all columns
split_and_save <- function(df, list_name, output_dir, train_ratio = 0.7, dev_ratio = 0.15) {
  
  # Remove duplicate rows based on the Sequence column
  df <- df %>% distinct(Sequence, .keep_all = TRUE)
  
  # Create output directory for this list
  dir_path <- file.path(output_dir, list_name)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Create subdirectory for full datasets
  full_dir_path <- file.path(dir_path, "all")
  dir.create(full_dir_path, recursive = TRUE, showWarnings = FALSE)
  
  # Stratified sampling
  train_indices <- createDataPartition(df$Label, p = train_ratio, list = FALSE)
  train <- df[train_indices, ]
  remaining <- df[-train_indices, ]
  
  dev_test_ratio <- dev_ratio / (1 - train_ratio)  # Adjust ratio for remaining data
  dev_indices <- createDataPartition(remaining$Label, p = dev_test_ratio, list = FALSE)
  dev <- remaining[dev_indices, ]
  test <- remaining[-dev_indices, ]
  
  # Write CSV files with only `Sequence` and `Label`
  write.csv(train[, c("Sequence", "Label")], file = file.path(dir_path, "train.csv"), row.names = FALSE)
  write.csv(dev[, c("Sequence", "Label")], file = file.path(dir_path, "dev.csv"), row.names = FALSE)
  write.csv(test[, c("Sequence", "Label")], file = file.path(dir_path, "test.csv"), row.names = FALSE)
  
  # Write CSV files with all columns
  write.csv(train, file = file.path(full_dir_path, "train_full.csv"), row.names = FALSE)
  write.csv(dev, file = file.path(full_dir_path, "dev_full.csv"), row.names = FALSE)
  write.csv(test, file = file.path(full_dir_path, "test_full.csv"), row.names = FALSE)
  
  cat("Data saved for:", list_name, "\n")
}


#  save
output_data_dnabert_dir <- paste0(save_pth_dna, "dnabert_2_data_nokmer/", bed_nm, "/unique_TSS/")
lapply(names(list_unique_TSS), function(name) {
  split_and_save(list_unique_TSS[[name]], name, output_data_dnabert_dir)
})

output_data_dnabert_dir <- paste0(save_pth_dna, "dnabert_2_data_nokmer/", bed_nm, "/unique_TSS_diff_50/")
lapply(names(list_unique_TSS_TSSdiff50), function(name) {
  split_and_save(list_unique_TSS_TSSdiff50[[name]], name, output_data_dnabert_dir)
})

cat("All files have been saved successfully.\n")

##--------------------------------------------------------------------------------------------------------
# data for shorter lengths compared to 10k TSS 2k
##--------------------------------------------------------------------------------------------------------

# saveRDS(list_unique_TSS_TSSdiff50, paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS_diff_50.rds"))
# saveRDS(list_unique_TSS, paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS.rds"))

# Function to extract strand information from a dataframe
extract_strand <- function(df) {
  df %>%
    mutate(
      strand = sapply(strsplit(Name, "\\|"), function(x) x[6]) # Extract 6th field
    )
}

bed_nm = "10k_TSS_2k_Correct" 
list_unique_TSS_TSSdiff50 = readRDS(paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS_diff_50.rds"))
list_unique_TSS =readRDS(paste0(save_pth_dna, "rds/", bed_nm, "/unique_TSS.rds"))

# Apply the function to each dataframe in the list
list_unique_TSS_TSSdiff50 <- lapply(list_unique_TSS_TSSdiff50, extract_strand)
list_unique_TSS <- lapply(list_unique_TSS, extract_strand)


# Range: start_modified = transcription_start_site - 10000
# end_modified   = transcription_start_site + 1999
# TSS Position: 10,001st base in the range
# 
# Range: start_modified = transcription_start_site - 1999
# end_modified   = transcription_start_site + 10000
# TSS Position: 1001st base in the reverse-complemented range

# Iterate over each dataframe in the list and apply the transformations
processed_list <- lapply(list_unique_TSS_TSSdiff50, function(data) {
  # Add a new column for the TSS position based on strand
  data <- data %>% 
    mutate(TSS_position = ifelse(strand == "+", 10001, 1001))
  
  # Apply each range cut and store in a list
  range_cuts <- list(
    "500_TSS_50" = c(-500, 50),
    "1k_TSS_200" = c(-1000, 200),
    "1k_TSS_1k" = c(-1000, 1000),
    "2k_TSS_1k" = c(-2000, 1000),
    "2k_TSS_500" = c(-2000, 500),
    "3k_TSS_1k" = c(-3000, 1000),
    "4k_TSS_1k" = c(-4000, 1000)
  )
  
  # Process each range
  range_results <- lapply(names(range_cuts), function(range_name) {
    range1 <- range_cuts[[range_name]]
    
    data <- data %>%
      mutate(
        start_idx = ifelse(strand == "+",
                           as.integer(TSS_position) + as.integer(range1[1]),
                           as.integer(TSS_position) - as.integer(range1[2])),
        end_idx = ifelse(strand == "+",
                         as.integer(TSS_position) + as.integer(range1[2]),
                         as.integer(TSS_position) - as.integer(range1[1]))
      ) %>%
      # Ensure indices are within bounds
      mutate(
        start_idx = pmax(1, start_idx),
        end_idx = pmin(nchar(Sequence), end_idx)
      )
    
    # Extract substring
    data <- data %>%
      mutate(
        Sequence = substr(Sequence, start_idx, end_idx)
      )
    
    return(data)
  })
  
  # Name each element in the range_results list
  names(range_results) <- names(range_cuts)
  
  # Return the list of processed ranges
  range_results
})


# SAVE

# Rename list function
rename_list <- function(old_names) {
  new_names <- old_names %>%
    tolower() %>% # Convert to lowercase
    gsub("^([^_]+).*_([^_]+)$", "\\1_\\2", .) %>% # Keep the first and last pieces separated by "_"
    gsub("trans", "", .) # Remove "trans" from the names
  return(new_names)
}
rename_list(names(processed_list)) -> names(processed_list)


# Create the main directories and subdirectories
base_dir <- paste0(save_pth_dna, "dnabert_2_data_nokmer", "/unique_TSS_diff_50/")
dir.create(base_dir, showWarnings = FALSE)

# Iterate over the names of `processed_list[[1]]` for main directories
main_dirs <- names(processed_list[[1]])

for (main_dir in main_dirs) {
  main_path <- file.path(base_dir, main_dir)
  dir.create(main_path, showWarnings = FALSE)
  
  # Iterate over `processed_list` names for subdirectories
  sub_dirs <- names(processed_list)
  for (sub_dir in sub_dirs) {
    sub_path <- file.path(main_path, sub_dir)
    dir.create(sub_path, showWarnings = FALSE)
    
    # Perform stratified sampling and save files if data exists
    if (!is.null(processed_list[[sub_dir]][[main_dir]])) {
      split_and_save(
        df = processed_list[[sub_dir]][[main_dir]],
        list_name = sub_dir,
        output_dir = main_path,
        train_ratio = 0.7,
        dev_ratio = 0.15
      )
    }
  }
}






