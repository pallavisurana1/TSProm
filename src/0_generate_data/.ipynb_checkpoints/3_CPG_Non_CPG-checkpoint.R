

# R4.4

# Load required libraries
pacman::p_load(data.table, dplyr, biomaRt, GenomicRanges, Biostrings, BSgenome, stringr)

# BiocManager::install(c(                                                
#   "BSgenome.Hsapiens.UCSC.hg38",
#   "BSgenome.Mmusculus.UCSC.mm39"
# ))

# pacman::p_load(BSgenome.Hsapiens.UCSC.hg38, BSgenome.Mmusculus.UCSC.mm39, BSgenome)

##--------------------------------------------------------------------------------------------------------
## Define paths
##--------------------------------------------------------------------------------------------------------

ft = "/Users/pallavisurana/Desktop/c/home/psurana"
ft = "/home/psurana"
save_pth_dna = paste0(ft, "/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/")
output_fasta_dir <- paste0(save_pth_dna, "fasta/")

# dir.create(output_fasta_dir, recursive = TRUE, showWarnings = FALSE)
# dir.create( paste0(save_pth_dna, "rds/"), recursive = TRUE, showWarnings = FALSE)

##--------------------------------------------------------------------------------------------------------
## Load Human & Mouse Data Separately
##--------------------------------------------------------------------------------------------------------
# 
# df_human = fread(file.path(ft, "/projects/Davuluri_lab/TransTEx_GTEX_Feb_2025/TransTEx_Normal_GTEx_gene_ENSEMBL113.csv"))
# df_mouse = readRDS(file.path(ft, "/projects/Davuluri_lab/Extramapper_2024/TransTEx/MouseTranstex/mouse_transtex_biomart.rds"))
# 
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
# 
# # Assign species labels
# df_human$species <- "human"
# df_mouse$species <- "mouse"
# 
# # Keep species-specific chromosomes
# df_human <- df_human %>% filter(chromosome_name %in% c(as.character(1:22), "X", "Y"))
# df_mouse <- df_mouse %>% filter(chromosome_name %in% c(as.character(1:19), "X", "Y"))
# 
# # Select required columns
# df_human1 <- df_human %>% dplyr::select(chromosome_name, transcript_start, transcript_end, transcription_start_site, strand, ensembl_gene_id, ensembl_transcript_id, tissue, category_final, species)
# df_mouse1 <- df_mouse %>% dplyr::select(chromosome_name, transcript_start, transcript_end, transcription_start_site, strand, transcript, tissue, category_final, species)
# 
# df_human1 <- df_human1 %>% distinct(transcription_start_site, .keep_all = TRUE) %>% dplyr::rename(transcript = ensembl_transcript_id)
# df_mouse1 <- df_mouse1 %>% distinct(transcription_start_site, .keep_all = TRUE)
# 
# dim(df_human1)
# dim(df_mouse1)
# 
# ##--------------------------------------------------------------------------------------------------------
# ## Compute BED Coordinates for Sequence Extraction
# ##--------------------------------------------------------------------------------------------------------
# 
# compute_bed_regions <- function(df) {
#   df %>%
#     mutate(
#       chromosome_name = paste0("chr", chromosome_name),
#       strand = ifelse(strand == 1, "+", "-"),
#       
#       `1k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 1000, transcription_start_site - 1000),
#       `1k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 1000),
#       
#       `2k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 2000, transcription_start_site - 1000),
#       `2k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 2000),
#       
#       `3k_TSS_1k_start` = ifelse(strand == "+", transcription_start_site - 3000, transcription_start_site - 1000),
#       `3k_TSS_1k_end`   = ifelse(strand == "+", transcription_start_site + 1000, transcription_start_site + 3000)
#     )
# }
# 
# df_human1 <- compute_bed_regions(df_human1)
# df_mouse1 <- compute_bed_regions(df_mouse1)
# 
# 
# ##--------------------------------------------------------------------------------------------------------
# ## Convert BED to FASTA Separately for Human and Mouse
# ##--------------------------------------------------------------------------------------------------------
# 
# # Define genomes
# genome_human <- BSgenome.Hsapiens.UCSC.hg38
# genome_mouse <- BSgenome.Mmusculus.UCSC.mm39
# 
# # Function to extract sequences
# extract_sequences <- function(df, genome) {
#   # List of columns containing start and end regions
#   region_types <- c("1k_TSS_1k", "2k_TSS_1k", "3k_TSS_1k")
#   
#   # Loop through each region type and extract sequences
#   for (region in region_types) {
#     start_col <- paste0(region, "_start")
#     end_col <- paste0(region, "_end")
#     
#     # Convert dataframe to GRanges
#     gr <- GRanges(
#       seqnames = df$chromosome_name,
#       ranges = IRanges(start = df[[start_col]], end = df[[end_col]]),
#       strand = df$strand
#     )
#     
#     # Extract sequences
#     sequences <- getSeq(genome, gr)
#     
#     # Add sequences to dataframe
#     df[[paste0(region, "_seq")]] <- as.character(sequences)
#   }
#   
#   return(df)
# }
# 
# df_mouse2 <- extract_sequences(df_mouse1, genome_mouse)
# df_human2 <- extract_sequences(df_human1, genome_human)
# 
# saveRDS(df_human2, paste0(save_pth_dna, "rds/", "human_seq.rds"))
# saveRDS(df_mouse2, paste0(save_pth_dna, "rds/", "mouse_seq.rds"))

##--------------------------------------------------------------------------------------------------------
# fetch fasta
##--------------------------------------------------------------------------------------------------------

df_human2 = readRDS(paste0(save_pth_dna, "rds/", "human_seq.rds"))
df_mouse2 = readRDS(paste0(save_pth_dna, "rds/", "mouse_seq.rds"))

##--------------------------------------------------------------------------------------------------------
## Merge Human & Mouse FASTA After Extraction and Perform Train-Dev-Test Split
##--------------------------------------------------------------------------------------------------------

# df_human2
# df_mouse2

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

##--------------------------------------------------------------------------------------------------------
# calculate cpg score
##--------------------------------------------------------------------------------------------------------
human_groups <- split(df_human2, df_human2$category_final)
mouse_groups <- split(df_mouse2, df_mouse2$category_final)

# Define valid tissue-category pairs
tissues_of_i = c("brain", "testis", "liver")
human_groups[["tsp"]] = human_groups[["tsp"]] %>% filter(tissue %in% tissues_of_i)
mouse_groups[["tsp"]] = mouse_groups[["tsp"]] %>% filter(tissue %in% tissues_of_i)


##--------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------


# Optimized function for generating subsequences
generate_subsequences_fast <- function(header, sequence, window_size = 201, step_size = 1) {
  seq_len <- nchar(sequence)
  starts <- seq(1, seq_len - window_size + 1, by = step_size)
  
  if (length(starts) == 0) return(NULL)  # Handle short sequences
  
  # Preallocate a data.table for speed
  dt <- data.table(
    Header = rep(header, length(starts)),
    Subseq_Num = seq_along(starts),
    Start = starts,
    End = starts + window_size - 1,
    Subsequence = substring(sequence, starts, starts + window_size - 1)
  )
  
  return(dt)
}

sequence_cols <- c("1k_TSS_1k_seq", "2k_TSS_1k_seq", "3k_TSS_1k_seq")


mouse_groups %>% names
# [1] "low"  "null" "tenh" "tsp"  "wide"
human_groups %>% names
# [1] "low"  "null" "tenh" "tsp"  "wide"

##-------------------------------------------
# Generate subsequences for all rows efficiently
##-------------------------------------------
mouse_subseq_results <- list()
human_subseq_results <- list()

# Function to process a group (shared by mouse and human)
process_groups <- function(group_list, species_label = "mouse") {
  result_list <- list()
  
  for (group_name in names(group_list)) {
    df <- group_list[[group_name]]
    
    for (seq_col in sequence_cols) {
      region_name <- gsub("_seq$", "", seq_col)  # e.g., "1k_TSS_1k"
      
      subseq_dt <- rbindlist(lapply(1:nrow(df), function(i) {
        sequence <- df[[seq_col]][i]
        if (is.na(sequence) || nchar(sequence) < 201) return(NULL)
        header <- paste(df$chromosome_name[i], df[[paste0(region_name, "_start")]][i], df[[paste0(region_name, "_end")]][i], sep = "_")
        generate_subsequences_fast(header, sequence)
      }), use.names = TRUE, fill = TRUE)
      
      # Add metadata
      subseq_dt[, Region := region_name]
      subseq_dt[, Group := group_name]
      subseq_dt[, Species := species_label]
      setcolorder(subseq_dt, c("Species", "Group", "Region", "Header", "Subseq_Num", "Start", "End", "Subsequence"))
      
      # Save to list
      result_list[[paste(species_label, group_name, region_name, sep = "_")]] <- subseq_dt
    }
  }
  
  return(result_list)
}

# Apply to both species
mouse_subseq_results <- process_groups(mouse_groups, species_label = "mouse")
human_subseq_results <- process_groups(human_groups, species_label = "human")

##--------------------------------------------------------------------------------------------------------

# Function to adjust coordinates
adjust_coordinates <- function(df) {
  df <- df %>%
    rowwise() %>%
    mutate(Chromosome = str_split(Header, "\\_")[[1]][1],
           Start = as.numeric(str_split(Header, "\\_")[[1]][2]) + Start,
           End = as.numeric(str_split(Header, "\\_")[[1]][2]) + End)
  return(df)
}

# Apply to mouse
mouse_subseq_results <- lapply(mouse_subseq_results, adjust_coordinates)

# Apply to human
human_subseq_results <- lapply(human_subseq_results, adjust_coordinates)


##--------------------------

# Function to calculate CpG percentage
calculate_cpg_percentage <- function(sequence) {
  sequence <- toupper(sequence)
  num_CpG <- str_count(sequence, "CG")
  window_length <- nchar(sequence)
  
  if (window_length > 1) {
    return((num_CpG / (window_length - 1)) * 100)
  } else {
    return(0)
  }
}

# Compute CpG percentage for each subsequence

# For mouse
mouse_subseq_results1 <- lapply(mouse_subseq_results, function(dt) {
  setDT(dt)
  dt[, CpG_Percentage := sapply(Subsequence, calculate_cpg_percentage)]
  return(dt)
})

# For human
human_subseq_results1 <- lapply(human_subseq_results, function(dt) {
  setDT(dt)
  dt[, CpG_Percentage := sapply(Subsequence, calculate_cpg_percentage)]
  return(dt)
})

##--------------------------
# Compute max CpG percentage and corresponding window for each sequence
##--------------------------

get_max_cpg_classification <- function(dt, list_name, cutoff = 6.5) {
  # Determine region window size and CpG region of interest
  if (grepl("1k_TSS_1k$", list_name)) {
    region_bp <- 2000
    region_start <- 501
    region_end <- 1500
  } else if (grepl("2k_TSS_1k$", list_name)) {
    region_bp <- 3000
    region_start <- 1501
    region_end <- 2500
  } else if (grepl("3k_TSS_1k$", list_name)) {
    region_bp <- 4000
    region_start <- 2501
    region_end <- 3500
  } else {
    stop("Unrecognized region pattern in list name: ", list_name)
  }
  
  # Determine which subsequence chunks fall within the region of interest
  chunk_starts <- seq(1, region_bp, by = 201)
  chunk_ends <- pmin(chunk_starts + 200, region_bp)
  which_chunks <- which(chunk_starts <= region_end & chunk_ends >= region_start)
  
  # Get the subsequence with the max CpG percentage per original region
  dt_max <- dt[, .SD[which.max(CpG_Percentage)], by = Header]
  
  # Classify based on CpG cutoff and location within the desired region
  dt_max[, Category := ifelse(CpG_Percentage >= cutoff & Subseq_Num %in% which_chunks,
                              "High CpG", "Low CpG")]
  
  # Return selected columns
  return(dt_max[, .(Chromosome, Start, End, Header, CpG_Percentage, Subsequence, Subseq_Num, Category)])
}


# Apply to mouse list
mouse_max_cpg_results <- mapply(
  get_max_cpg_classification,
  dt = mouse_subseq_results1,
  list_name = names(mouse_subseq_results1),
  SIMPLIFY = FALSE
)

# Apply to human list
human_max_cpg_results <- mapply(
  get_max_cpg_classification,
  dt = human_subseq_results1,
  list_name = names(human_subseq_results1),
  SIMPLIFY = FALSE
)

save_dir_rds="/home/psurana/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/rds/cpg_n_cpg"
saveRDS(mouse_max_cpg_results, file.path(save_dir_rds, "mouse_raw1.rds"))
saveRDS(human_max_cpg_results, file.path(save_dir_rds, "human_raw1.rds"))

##--------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------
##--------------------------------------------------------------------------------------------------------


