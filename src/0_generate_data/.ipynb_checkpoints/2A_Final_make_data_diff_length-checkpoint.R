

# R411
.libPaths("/vast/projects/rdavuluri-group/psurana/conda_environments/r451/lib/R/library")

# Load required libraries
pacman::p_load(data.table, dplyr, biomaRt, GenomicRanges, Biostrings, BSgenome, 
               BSgenome.Hsapiens.UCSC.hg38, BSgenome.Mmusculus.UCSC.mm39, caret)
##--------------------------------------------------------------------------------------------------------
## Define paths
##--------------------------------------------------------------------------------------------------------

ft = "/Users/pallavisurana/Desktop/c/home/psurana"
ft = "/vast/projects/rdavuluri-group/psurana/"
save_pth_dna = paste0(ft, "/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/")
output_fasta_dir <- paste0(save_pth_dna, "fasta/")

dir.create(output_fasta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create( paste0(save_pth_dna, "rds/"), recursive = TRUE, showWarnings = FALSE)

##--------------------------------------------------------------------------------------------------------
## Load Human & Mouse Data Separately
##--------------------------------------------------------------------------------------------------------

df_human = fread(file.path(ft, "/projects/Davuluri_lab/TransTEx_GTEX_Feb_2025/TransTEx_Normal_GTEx_gene_ENSEMBL113.csv"))
df_mouse = readRDS(file.path(ft, "/projects/Davuluri_lab/Extramapper_2024/TransTEx/MouseTranstex/mouse_transtex_biomart.rds"))

# ##--------------------------------------------------------------------------------------------------------
# ## Retrieve Chromosome Locations for Mouse Transcripts (Missing Data)
# ##--------------------------------------------------------------------------------------------------------
# 
# df_mouse = readRDS(file.path(ft, "/projects/Davuluri_lab/Extramapper_2024/TransTEx/MouseTranstex/mouse_transtex_modified.rds"))
# df_mouse <- df_mouse %>% dplyr::select(transcript, `.id`, tissue) %>% dplyr::rename(category_final = `.id`)
# 
# # Connect to Ensembl BioMart
# ensembl_mouse <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
#                          dataset = "mmusculus_gene_ensembl", 
#                          host = "https://www.ensembl.org")
# 
# # Query for missing chromosome locations based on ENSEMBL Transcript IDs
# mouse_location_data <- getBM(attributes = c("ensembl_transcript_id", "chromosome_name", 
#                                             "transcript_start", "transcript_end", 
#                                             "strand", "transcription_start_site"),
#                              filters = "ensembl_transcript_id",
#                              values = (df_mouse$transcript),
#                              mart = ensembl_mouse)
# 
# # Merge with df_mouse
# df_mouse1 <- merge(df_mouse, mouse_location_data, by.x = "transcript", by.y = "ensembl_transcript_id")
# saveRDS(df_mouse1, file.path(ft, "/projects/Davuluri_lab/Extramapper_2024/TransTEx/MouseTranstex/mouse_transtex_biomart.rds"))
# 
# ##--------------------------------------------------------------------------------------------------------

# Assign species labels
df_human$species <- "human"
df_mouse$species <- "mouse"

# Keep species-specific chromosomes
df_human <- df_human %>% filter(chromosome_name %in% c(as.character(1:22), "X", "Y"))
df_mouse <- df_mouse %>% filter(chromosome_name %in% c(as.character(1:19), "X", "Y"))

# Select required columns
df_human1 <- df_human %>% dplyr::select(chromosome_name, transcript_start, transcript_end, transcription_start_site, strand, ensembl_gene_id, ensembl_transcript_id, tissue, category_final, species)
df_mouse1 <- df_mouse %>% dplyr::select(chromosome_name, transcript_start, transcript_end, transcription_start_site, strand, transcript, tissue, category_final, species)

df_human1 <- df_human1 %>% distinct(transcription_start_site, .keep_all = TRUE) %>% dplyr::rename(transcript = ensembl_transcript_id)
df_mouse1 <- df_mouse1 %>% distinct(transcription_start_site, .keep_all = TRUE)

dim(df_human1)
dim(df_mouse1)

# just check tsp ones for some summary
df_mouse_tsp = df_mouse1 %>% filter(category_final == "tsp"); df_mouse_tsp$tissue %>% table %>% as.data.frame
df_human_tsp = df_human1 %>% filter(category_final == "TspTrans"); df_human_tsp$tissue %>% table %>% as.data.frame

##--------------------------------------------------------------------------------------------------------
## Compute BED Coordinates for Sequence Extraction
##--------------------------------------------------------------------------------------------------------

compute_bed_regions <- function(df) {
  df %>%
    mutate(
      chromosome_name = paste0("chr", chromosome_name),
      strand = ifelse(strand == 1, "+", "-"),
      
      # Define TSS-centered regions
      `500_TSS_50_start` = ifelse(strand == "+", transcription_start_site - 500, transcription_start_site - 50),
      `500_TSS_50_end`   = ifelse(strand == "+", transcription_start_site + 50, transcription_start_site + 500),
      
      `1k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 1000, transcription_start_site - 1000),
      `1k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 1000),
      
      `2k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 2000, transcription_start_site - 1000),
      `2k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 2000),
      
      `3k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 3000, transcription_start_site - 1000),
      `3k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 3000)
    )
}


df_human1 <- compute_bed_regions(df_human1)
df_mouse1 <- compute_bed_regions(df_mouse1)


##--------------------------------------------------------------------------------------------------------
## Convert BED to FASTA Separately for Human and Mouse
##--------------------------------------------------------------------------------------------------------

# Define genomes
genome_human <- BSgenome.Hsapiens.UCSC.hg38
genome_mouse <- BSgenome.Mmusculus.UCSC.mm39

# Function to extract sequences
extract_sequences <- function(df, genome) {
  # List of columns containing start and end regions
  region_types <- c("500_TSS_50", "1k_TSS_1k", "2k_TSS_1k", "3k_TSS_1k")
  
  # Loop through each region type and extract sequences
  for (region in region_types) {
    start_col <- paste0(region, "_start")
    end_col <- paste0(region, "_end")
    
    # Convert dataframe to GRanges
    gr <- GRanges(
      seqnames = df$chromosome_name,
      ranges = IRanges(start = df[[start_col]], end = df[[end_col]]),
      strand = df$strand
    )
    
    # Extract sequences
    sequences <- getSeq(genome, gr)
    
    # Add sequences to dataframe
    df[[paste0(region, "_seq")]] <- as.character(sequences)
  }
  
  return(df)
}

df_mouse2 <- extract_sequences(df_mouse1, genome_mouse)
df_human2 <- extract_sequences(df_human1, genome_human)

saveRDS(df_human2, paste0(save_pth_dna, "rds/", "human_seq.rds"))
saveRDS(df_mouse2, paste0(save_pth_dna, "rds/", "mouse_seq.rds"))

##--------------------------------------------------------------------------------------------------------
## Merge Human & Mouse FASTA After Extraction and Perform Train-Dev-Test Split
##--------------------------------------------------------------------------------------------------------

# df_human2
# df_mouse2

df_human2 = readRDS(paste0(save_pth_dna, "rds/", "human_seq.rds"))
df_mouse2 = readRDS(paste0(save_pth_dna, "rds/", "mouse_seq.rds"))
        
# Rename category_final values
df_human2 <- df_human2 %>%
  mutate(category_final = case_when(
    category_final == "TspTrans"  ~ "tsp",
    category_final == "NullTrans" ~ "null",
    category_final == "TenhTrans" ~ "tenh",
    category_final == "WideTrans" ~ "wide",
    category_final == "LowTrans"  ~ "low",
    TRUE ~ category_final  # Keep unchanged values
  ))  
df_human2$tissue=tolower(df_human2$tissue)


# Define tissues to keep within TSP group
selected_tissues <- c("brain", "testis", "liver")
selected_tissues <- c("muscle", "spleen")


# Function to match specific tissues within TSP to other categories
match_tsp_tissues <- function(df_human, df_mouse) {
  # Split by category_final for both species
  human_groups <- split(df_human, df_human$category_final)
  mouse_groups <- split(df_mouse, df_mouse$category_final)
  
  # Initialize an empty list to store category pairings
  matched_data <- list()
  
  # Define valid tissue-category pairs
  # valid_pairs <- list(
  #   "brain" = "tsp_brain",
  #   "testis" = "tsp_testis",
  #   "liver" = "tsp_liver"
  # )
  valid_pairs <- list(
    "muscle" = "tsp_muscle",
    "spleen" = "tsp_spleen"
    )
  
  # Categories to pair with TSP tissues
  other_categories <- c("null", "low", "tenh", "wide")
  
  # Iterate through each selected TSP tissue and pair it with other groups
  for (tissue in names(valid_pairs)) {
    tsp_group_name <- valid_pairs[[tissue]]  # Get correct TSP group name
    
    for (category in other_categories) {
      pair_name <- paste0(tsp_group_name, "_", category)  # Naming convention
      
      # Initialize an empty dataframe for this pair
      matched_df <- data.frame()
      
      # Match in human data
      if ("tsp" %in% names(human_groups) && category %in% names(human_groups)) {
        tsp_df <- human_groups[["tsp"]] %>%
          filter(tissue == tissue) %>%  # Keep only the relevant TSP tissue
          mutate(Label = 1)  # Ensure Label = 1
        
        other_df <- human_groups[[category]] %>%
          mutate(Label = 0)  # Ensure Label = 0
        
        matched_df <- bind_rows(matched_df, tsp_df, other_df)
      }
      
      # Match in mouse data
      if ("tsp" %in% names(mouse_groups) && category %in% names(mouse_groups)) {
        tsp_df <- mouse_groups[["tsp"]] %>%
          filter(tissue == tissue) %>%  # Keep only the relevant TSP tissue
          mutate(Label = 1)  # Ensure Label = 1
        
        other_df <- mouse_groups[[category]] %>%
          mutate(Label = 0)  # Ensure Label = 0
        
        matched_df <- bind_rows(matched_df, tsp_df, other_df)
      }
      
      # Store the matched dataframe for this TSP-category pair
      if (nrow(matched_df) > 0) {
        matched_data[[pair_name]] <- matched_df
      }
    }
  }
  
  return(matched_data)
}

matched_tsp_tissue_groups = match_tsp_tissues(df_human2, df_mouse2)

##------------------------------

# Define function to filter elements of the list
filter_tissue_data <- function(df, df_name) {
  # Determine tissue type based on list element name
  if (grepl("spleen", df_name)) {
    tissue_type <- "spleen"
  } else if (grepl("muscle", df_name)) {
    tissue_type <- "muscle"
  } else {
    return(df)  # Return unchanged if not in defined categories
  }
  
  # Apply filtering
  df %>%
    filter((category_final != "tsp") | (category_final == "tsp" & tissue == tissue_type))
}

# Apply function to each element in the list with corresponding names
matched_tsp_tissue_groups <- mapply(filter_tissue_data, matched_tsp_tissue_groups, names(matched_tsp_tissue_groups), SIMPLIFY = FALSE)


##--------------------------------------------------------------------------------------------------------
# Function to compute sequence length and categorize sequences by length within each tissue-category pair
##--------------------------------------------------------------------------------------------------------

# Function to compute sequence lengths for all `_seq` columns
calculate_sequence_lengths <- function(df) {
  # Identify columns that end with `_seq`
  seq_cols <- grep("_seq$", colnames(df), value = TRUE)
  
  # Compute sequence lengths for each sequence column
  for (seq_col in seq_cols) {
    length_col <- gsub("_seq$", "_length", seq_col)  # Create new column name
    df[[length_col]] <- nchar(df[[seq_col]])  # Compute length
  }
  
  return(df)
}

# Apply function to all TSP-tissue-category pairs
for (pair_name in names(matched_tsp_tissue_groups)) {
  matched_tsp_tissue_groups[[pair_name]] <- calculate_sequence_lengths(matched_tsp_tissue_groups[[pair_name]])
}
# saveRDS(matched_tsp_tissue_groups, paste0(save_pth_dna, "rds/", "2_matched_tsp_tissue_groups.rds"))
saveRDS(matched_tsp_tissue_groups, paste0(save_pth_dna, "rds/", "2_matched_tsp_spleen_tissue_groups_muscle_spleen.rds"))

##--------------------------------------------------------------------------------------------------------
# save fastas for mmseqs2 and then split
##--------------------------------------------------------------------------------------------------------
library(readr)
library(stringr)

# Combine each group/length into a single CSV with two columns: header, sequence

output_combined_dir <- paste0(save_pth_dna, "b4Split_combinedCSV/")
dir.create(output_combined_dir, recursive = TRUE, showWarnings = FALSE)

for (group in names(matched_tsp_tissue_groups)) {
  df <- matched_tsp_tissue_groups[[group]]
  if (!all(df$category_final == "tsp")) {
    df$tissue <- "other"
  }
  seq_cols <- grep("_seq$", names(df), value = TRUE)
  if (length(seq_cols) == 0) next  # Skip if no sequences
  
  for (seq_col in seq_cols) {
    # Build the header for every row
    headers <- paste(
      df$chromosome_name, df$strand, df$tissue, df$category_final,
      df$species, df$Label, df$transcription_start_site, df$transcript,
      sep = "|"
    )
    # Build output df
    out_df <- data.frame(
      header = headers,
      sequence = df[[seq_col]]
    )
    # Remove NAs or empty sequences
    out_df <- out_df[!is.na(out_df$sequence) & out_df$sequence != "", ]
    # Write CSV, e.g. group_500_TSS_50_seq.csv
    out_file <- file.path(output_combined_dir, paste0(group, "_", seq_col, ".csv"))
    write.csv(out_df, out_file, row.names = FALSE, quote = TRUE)
    cat("Wrote:", out_file, "\n")
  }
}

# ##--------------------------------------------------------------------------------------------------------
# # 1:1 split then make test train and dev from it
# ##--------------------------------------------------------------------------------------------------------
# 
# # Define sequence length groups
# length_groups <- c(550, 2000, 3000, 4000)
# 
# # Initialize list to store balanced data split by length
# balanced_split_by_length <- list()
# 
# # Function to split data based on length and sequence columns, while maintaining balance
# split_and_balance_data <- function(df, category_name) {
#   # Identify sequence and length columns dynamically
#   seq_cols <- grep("_seq$", colnames(df), value = TRUE)
#   length_cols <- grep("_length$", colnames(df), value = TRUE)
#   
#   # Ensure mapping between sequences and their respective length columns
#   length_map <- setNames(length_cols, seq_cols)
#   
#   # Create a list for storing length-based splits
#   category_splits <- list()
#   
#   for (seq_col in seq_cols) {
#     length_col <- length_map[[seq_col]]  # Get corresponding length column
#     
#     for (len in length_groups) {
#       subset_df <- df %>%
#         dplyr::filter(.data[[length_col]] == len+1) %>%
#         dplyr::select(chromosome_name, strand, tissue, category_final, species, Label, seq_col, length_col, transcription_start_site, transcript)
#       
#       # Ensure there are both Label = 1 and Label = 0 samples
#       if (nrow(subset_df) > 0) {
#         count_label_1 <- nrow(subset_df %>% filter(Label == 1))
#         count_label_0 <- nrow(subset_df %>% filter(Label == 0))
#         
#         if (count_label_1 > 0 & count_label_0 > 0) {
#           min_count <- min(count_label_1, count_label_0)  # Balance dataset
#           
#           balanced_subset <- bind_rows(
#             subset_df %>% filter(Label == 1) %>% sample_n(min_count, replace = FALSE),
#             subset_df %>% filter(Label == 0) %>% sample_n(min_count, replace = FALSE)
#           )
#           
#           split_name <- paste0(category_name, "_", len)
#           category_splits[[split_name]] <- balanced_subset
#         }
#       }
#     }
#   }
#   
#   return(category_splits)
# }
# 
# # Apply function to all matched TSP tissue groups
# for (category_name in names(matched_tsp_tissue_groups)) {
#   balanced_split_by_length[[category_name]] <- split_and_balance_data(matched_tsp_tissue_groups[[category_name]], category_name)
# }
# 
# 
# ##--------------------------------------------------------------------------------------------------------
# 
# # Function to split data into train, dev, and test based on labels
# split_data <- function(df) {
#   if (!("Label" %in% colnames(df))) return(df) # Ensure 'label' column exists
#   
#   # Rename the column ending with "_seq" to "Sequence"
#   df <- df %>%
#     dplyr::rename(Sequence = dplyr::select(., ends_with("_seq")) %>% names())
#   
#   # Select required columns and create a new column by concatenating specified ones
#   df <- df %>%
#     dplyr::select(Sequence, Label, chromosome_name, strand, tissue, category_final, species, transcription_start_site, transcript) %>%
#     dplyr::mutate(concatenated_column = paste(chromosome_name, strand, tissue, category_final, species, Label, transcription_start_site, transcript, sep = "|"))
#   
#   df_1 <- df %>% filter(Label == 1)
#   df_0 <- df %>% filter(Label == 0)
#   
#   train_1 <- df_1[1:floor(0.8 * nrow(df_1)), ]
#   dev_1 <- df_1[(floor(0.8 * nrow(df_1)) + 1):floor(0.9 * nrow(df_1)), ]
#   test_1 <- df_1[(floor(0.9 * nrow(df_1)) + 1):nrow(df_1), ]
#   
#   train_0 <- df_0[1:floor(0.8 * nrow(df_0)), ]
#   dev_0 <- df_0[(floor(0.8 * nrow(df_0)) + 1):floor(0.9 * nrow(df_0)), ]
#   test_0 <- df_0[(floor(0.9 * nrow(df_0)) + 1):nrow(df_0), ]
#   
#   list(
#     train = bind_rows(train_1, train_0),
#     dev = bind_rows(dev_1, dev_0),
#     test = bind_rows(test_1, test_0)
#   )
# }
# 
# # Apply split function to each sublist in the list of lists
# split_results <- lapply(balanced_split_by_length, function(sublist) {
#   lapply(sublist, split_data)
# })
# 
# # Define output directory
# output_split_dir <- paste0(save_pth_dna, "split/")
# output_splitPheno_dir <- paste0(save_pth_dna, "split/pheno/")
# 
# # Function to extract last part of names
# extract_length <- function(name) {
#   sub(".*_", "", name)  # Extracts the last numeric portion after the last underscore
# }
# 
# # Save split data
# for (group in names(split_results)) {
#   for (len in names(split_results[[group]])) {
#     
#     len_extracted <- extract_length(len)  # Extract numeric part
#     
#     dir.create(file.path(output_split_dir, group, len_extracted), recursive = TRUE, showWarnings = FALSE)
#     dir.create(file.path(output_splitPheno_dir, group, len_extracted), recursive = TRUE, showWarnings = FALSE)
#     
#     write.csv(split_results[[group]][[len]]$train[ , c("Sequence", "Label")], file = paste0(output_split_dir, group, "/", len_extracted, "/train.csv"), row.names = FALSE)
#     write.csv(split_results[[group]][[len]]$dev[ , c("Sequence", "Label")], file = paste0(output_split_dir, group, "/", len_extracted, "/dev.csv"), row.names = FALSE)
#     write.csv(split_results[[group]][[len]]$test[ , c("Sequence", "Label")], file = paste0(output_split_dir, group, "/", len_extracted, "/test.csv"), row.names = FALSE)
#     
#     write.csv(split_results[[group]][[len]]$train[ , c("concatenated_column")], file = paste0(output_splitPheno_dir, group, "/", len_extracted, "/train.csv"), row.names = FALSE)
#     write.csv(split_results[[group]][[len]]$dev[ , c("concatenated_column")], file = paste0(output_splitPheno_dir, group, "/", len_extracted, "/dev.csv"), row.names = FALSE)
#     write.csv(split_results[[group]][[len]]$test[ , c("concatenated_column")], file = paste0(output_splitPheno_dir, group, "/", len_extracted, "/test.csv"), row.names = FALSE)
#   }
# }









###------------------------------------------------------------------------------------------------------------
## July 4, 2024
###------------------------------------------------------------------------------------------------------------
# make the 4 new types of data for fine tuning
# TSPall vs rest (Human only)
# TSPall vs Non-Prom (Human only)
# TSP site vs rest
# TSP site vs Non-Prom
###------------------------------------------------------------------------------------------------------------

df_human2 = readRDS(paste0(save_pth_dna, "rds/", "human_seq.rds"))
df_mouse2 = readRDS(paste0(save_pth_dna, "rds/", "mouse_seq.rds"))

# Rename category_final values
df_human2 <- df_human2 %>%
  mutate(category_final = case_when(
    category_final == "TspTrans"  ~ "tsp",
    category_final == "NullTrans" ~ "null",
    category_final == "TenhTrans" ~ "tenh",
    category_final == "WideTrans" ~ "wide",
    category_final == "LowTrans"  ~ "low",
    TRUE ~ category_final  # Keep unchanged values
  ))  
df_human2$tissue=tolower(df_human2$tissue)
df_human2$ensembl_gene_id = NULL

transtex = rbind(df_human2, df_mouse2)
transtex$category_final %>% table()
transtex$header <- paste(
  transtex$chromosome_name,
  transtex$transcription_start_site,
  transtex$strand,
  transtex$transcript,
  transtex$tissue,
  transtex$species,
  sep = "_"
)
transtex_tsp = transtex %>% filter(category_final == "tsp")
transtex_rest = transtex %>% filter(category_final != "tsp"); transtex_rest$Label = 0
transtex_rest$header <- paste(
  transtex_rest$chromosome_name,
  transtex_rest$transcription_start_site,
  transtex_rest$strand,
  transtex_rest$transcript,
  transtex_rest$species,
  sep = "_"
)
non_promoter_df = readRDS(paste0(save_pth_dna, "rds/hg38_nonPromoter_randomSeqs.rds"))
non_promoter_df$header = "non_prom_human"

save_pth = paste0(save_pth_dna, "july4_data_vs_NonPROM/")

###-------------------------
# TSp all vs non prom human
colnames_to_select = intersect(transtex_tsp %>% colnames, non_promoter_df %>% colnames)
non_promoter_df1 = non_promoter_df %>% select(all_of(colnames_to_select))
non_promoter_df1$Label = 0
transtex_tsp1 = transtex_tsp %>% select(all_of(colnames_to_select))
transtex_tsp1$Label = 1

tspAll_nonPromHu = rbind(transtex_tsp1, non_promoter_df1)
tspAll_nonPromHu$Label %>% table

write.csv(tspAll_nonPromHu, paste0(save_pth, "tspAll_nonPromHu.csv"))

###-------------------------
# TSp site vs non prom human
transtex_tsp1_brain = transtex_tsp %>% filter(tissue == "brain") %>% select(all_of(colnames_to_select))
transtex_tsp1_testis= transtex_tsp %>% filter(tissue == "testis")%>% select(all_of(colnames_to_select))
transtex_tsp1_liver= transtex_tsp %>% filter(tissue == "liver")%>% select(all_of(colnames_to_select))
transtex_tsp1_spleen= transtex_tsp %>% filter(tissue == "spleen")%>% select(all_of(colnames_to_select))
transtex_tsp1_muscle= transtex_tsp %>% filter(tissue == "muscle")%>% select(all_of(colnames_to_select))

transtex_tsp1_brain$Label = transtex_tsp1_liver$Label = transtex_tsp1_testis$Label = 1
transtex_tsp1_spleen$Label = transtex_tsp1_muscle$Label = 1

tspbrain_nonPromHu = rbind(non_promoter_df1, transtex_tsp1_brain)
tsptestis_nonPromHu = rbind(non_promoter_df1, transtex_tsp1_testis)
tspliver_nonPromHu = rbind(non_promoter_df1, transtex_tsp1_liver)
tspspleen_nonPromHu = rbind(non_promoter_df1, transtex_tsp1_spleen)
tspmuscle_nonPromHu = rbind(non_promoter_df1, transtex_tsp1_muscle)

write.csv(tspbrain_nonPromHu, paste0(save_pth, "tspbrain_nonPromHu.csv"))
write.csv(tsptestis_nonPromHu, paste0(save_pth, "tsptestis_nonPromHu.csv"))
write.csv(tspliver_nonPromHu, paste0(save_pth, "tspliver_nonPromHu.csv"))
write.csv(tspspleen_nonPromHu, paste0(save_pth, "tspspleen_nonPromHu.csv"))
write.csv(tspmuscle_nonPromHu, paste0(save_pth, "tspmuscle_nonPromHu.csv"))

###-------------------------
# TSp all vs Non TSp all
colnames_to_select = intersect(transtex_rest %>% colnames, transtex_tsp1 %>% colnames)
colnames_to_select

transtex_rest1 = transtex_rest %>% select(all_of(colnames_to_select))
tspAll_restTranstex = rbind(transtex_tsp1, transtex_rest1)

write.csv(tspAll_restTranstex, paste0(save_pth, "tspAll_restTranstex.csv"))

###-------------------------
# TSp site vs Non TSp all

tspbrain_restTranstex = rbind(transtex_rest1, transtex_tsp1_brain)
tsptestis_restTranstex = rbind(transtex_rest1, transtex_tsp1_testis)
tspliver_restTranstex = rbind(transtex_rest1, transtex_tsp1_liver)

write.csv(tspbrain_restTranstex, paste0(save_pth, "tspbrain_restTranstex.csv"))
write.csv(tsptestis_restTranstex, paste0(save_pth, "tsptestis_restTranstex.csv"))
write.csv(tspliver_restTranstex, paste0(save_pth, "tspliver_restTranstex.csv"))

###-------------------------###-------------------------###-------------------------
library(Biostrings)
library(readr)    
library(dplyr)
library(stringr)
library(tools)      

csv_dir   <- "/vast/projects/rdavuluri-group/psurana//projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/b4Split_combinedCSV/move_ramana_server"
csv_files <- list.files(csv_dir, pattern = "\\.csv$", full.names = TRUE)

csv_dir = save_pth
csv_files <- list.files(
  csv_dir,
  pattern = "(spleen|muscle).*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

write_fasta <- function(df, seq_col, out_fasta) {
  stopifnot(all(c("header", seq_col) %in% names(df)))
  dna <- DNAStringSet(df[[seq_col]])
  # names(dna) <- stringr::str_c(df$header, df$Label, sep = "|")
  names(dna) <- stringr::str_c(df$header)
  writeXStringSet(dna, filepath = out_fasta, format = "fasta")
}

for (csv in csv_files) {
  df <- data.table::fread(csv, data.table = FALSE, showProgress = FALSE)
  base <- file_path_sans_ext(basename(csv))
  out2k <- file.path(dirname(csv), paste0(base, "_2k.fasta"))
  out3k <- file.path(dirname(csv), paste0(base, "_3k.fasta"))
  out4k <- file.path(dirname(csv), paste0(base, "_4k.fasta"))
  
  write_fasta(df, "1k_TSS_1k_seq", out2k)
  write_fasta(df, "2k_TSS_1k_seq", out3k)
  write_fasta(df, "3k_TSS_1k_seq", out4k)
  
  message("wrote FASTAs for ", base)
}

 
# ################################################################################
# # get non prom data
# ## Random “background” sequences that sit ≥2 kb away from any TSS
# ## ────────────────────────────────────────────────────────────────────────────
# ## • Three region lengths: 2 000 bp, 3 000 bp, 4 000 bp
# ## • Their sequences are appended to your existing data frame           
# ##   as columns that match your original naming scheme:
# ##      – 1k_TSS_1k_seq  (for 2 000 bp random regions)
# ##      – 2k_TSS_1k_seq  (for 3 000 bp random regions)
# ##      – 3k_TSS_1k_seq  (for 4 000 bp random regions)
# ################################################################################
# 
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(GenomicFeatures)
# library(GenomicRanges)
# library(regioneR)
# library(dplyr)
# 
# genome <- BSgenome.Hsapiens.UCSC.hg38   # reference sequence object
# txdb   <- TxDb.Hsapiens.UCSC.hg38.knownGene
# 
# # ── 1.  Make a ±1 kb mask around every TSS ───────────────────────────────────
# tss       <- promoters(genes(txdb), upstream = 0, downstream = 1)   # 1-bp at TSS
# tss_mask  <- resize(tss, width = 2001, fix = "center")              # ±1 000 bp
# 
# # ── 2.  Decide how many random regions you want ──────────────────────────────
# n_regions <- 40000      # ← change this to any number you need
# set.seed(20250704)       # reproducible sampling
# # Helper function: generate random regions of fixed length away from TSS
# randSeqs <- function(len, n = n_regions) {
#   # Create dummy regions to be randomized
#   dummy <- GRanges(seqnames = sample(seqnames(genome), n, replace = TRUE),
#                    ranges = IRanges(start = rep(1, n), width = len))
# 
#   randomized <- randomizeRegions(A = dummy,
#                                  genome = genome,
#                                  mask = tss_mask,
#                                  allow.overlaps = FALSE,
#                                  per.chromosome = FALSE)
# 
#   seqs <- getSeq(genome, randomized)
#   as.character(seqs)
# }
# 
# # Generate sequences for 2k, 3k, 4k lengths
# seq_2k <- randSeqs(2000)
# seq_3k <- randSeqs(3000)
# seq_4k <- randSeqs(4000)
# 
# # ── 4.  Assemble one tidy data-frame ─────────────────────────────────────────
# non_promoter_df <- tibble(
#   id               = seq_len(n_regions),
#   `1k_TSS_1k_seq`  = seq_2k,
#   `2k_TSS_1k_seq`  = seq_3k,
#   `3k_TSS_1k_seq`  = seq_4k
# )
# # ── 5.  Inspect or save ──────────────────────────────────────────────────────
# saveRDS(non_promoter_df, paste0(save_pth_dna, "rds/hg38_nonPromoter_randomSeqs.rds"))

## ────────────────────────────────────────────────────── ──────────────────────────────────────────────────────
## run ~/projects/CODES_master/Davuluri_lab/DNABERT_TSpProm/Fine_Tune_dna2/human_mouse/2B_mmseqs.sh
## easy-cluster by default clusters the entries of a FASTA/FASTQ file using a cascaded clustering algorithm
## ────────────────────────────────────────────────────── ──────────────────────────────────────────────────────




