
pacman::p_load(data.table, dplyr, biomaRt, GenomicRanges, Biostrings, BSgenome, stringr)

##--------------------------------------------------------------------------------------------------------
## Define paths
##--------------------------------------------------------------------------------------------------------

ft = "/Users/pallavisurana/Desktop/c/home/psurana"
ft = "/home/psurana"
save_pth_dna = paste0(ft, "/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/")
output_fasta_dir <- paste0(save_pth_dna, "fasta/")

##--------------------------------------------------------------------------------------------------------
## Load Human & Mouse Data Separately; make cpg noncpg groups
##--------------------------------------------------------------------------------------------------------

# entire seq
df_human2 = readRDS(paste0(save_pth_dna, "rds/", "human_seq.rds"))
df_mouse2 = readRDS(paste0(save_pth_dna, "rds/", "mouse_seq.rds"))

paste(
  chromosome_name,
  transcript_start,
  transcript_end,
  transcription_start_site,
  strand,
  transcript,
  category_final,
  species,
  tissue,
  sep = "_" ))

# cpg_percent info
save_dir_rds="/home/psurana/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/rds/cpg_n_cpg"
mouse_max_cpg_results=readRDS(file.path(save_dir_rds, "mouse_raw.rds"))
human_max_cpg_results=readRDS(file.path(save_dir_rds, "human_raw.rds"))

##--------------------------------------------------------------------------------------------------------

# Get common suffixes (remove "human_" and "mouse_" from names)
common_suffixes <- gsub("^human_", "", names(human_max_cpg_results))

# Create merged list
merged_cpg_results <- lapply(common_suffixes, function(suffix) {
  human_data <- human_max_cpg_results[[paste0("human_", suffix)]]
  mouse_data <- mouse_max_cpg_results[[paste0("mouse_", suffix)]]
  rbind(human_data, mouse_data)
})

# Assign proper names
names(merged_cpg_results) <- common_suffixes


##--------------------------------------------------------------------------------------------------------
# Function to match specific tissues within TSP to other categories
match_tsp_tissues <- function(df_human, df_mouse) {
  # Split by category_final for both species
  human_groups <- split(df_human, df_human$category_final)
  mouse_groups <- split(df_mouse, df_mouse$category_final)
  
  # Initialize an empty list to store category pairings
  matched_data <- list()
  
  # Define valid tissue-category pairs
  valid_pairs <- list(
    "brain" = "tsp_brain",
    "testis" = "tsp_testis",
    "liver" = "tsp_liver"
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
  if (grepl("brain", df_name)) {
    tissue_type <- "brain"
  } else if (grepl("testis", df_name)) {
    tissue_type <- "testis"
  } else if (grepl("liver", df_name)) {
    tissue_type <- "liver"
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


##--------------------------------------------------------------------------------------------------------
# make some plots
##--------------------------------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)  # for arranging multiple plots

# 1) Your region suffixes
regions <- c("1k_TSS_1k", "2k_TSS_1k", "3k_TSS_1k")

# 2) For each region, bind the four category‐specific data.frames
combined_by_region <- lapply(regions, function(reg) {
  # find all names in merged_cpg_results ending with _<reg>
  pattern <- paste0("_", reg, "$")
  keys    <- grep(pattern, names(merged_cpg_results), value = TRUE)
  
  # bind them, tagging category
  df_list <- lapply(keys, function(k) {
    cat   <- sub(pattern, "", k)      # e.g. "low", "null", ...
    df    <- matched_tsp_tissue_groups[[k]]
    df$Category <- factor(cat, levels = c("low", "null", "wide", "tenh", "tsp"))
    df
  })
  combined <- do.call(rbind, df_list)
  combined$Region <- reg
  combined
})

names(combined_by_region) <- regions

# 3) Build one density plot per region
plots <- lapply(regions, function(reg) {
  df <- combined_by_region[[reg]]
  
  ggplot(df, aes(x = CpG_Percentage, color = Category, fill = Category)) +
    geom_density(alpha = 0.1, size = 1) +
    scale_color_manual(values = c(
      low  = "#1b9e77",
      null = "#d95f02",
      wide = "#7570b3",
      tenh = "#e7298a",
      tsp = "darkred"
    )) +
    labs(
      title = paste("Region", reg),
      x     = "CpG Percentage",
      y     = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title     = element_text(hjust = 0.5),
      legend.position= "bottom"
    )
})

# 4) Arrange the three plots side‐by‐side
grid.arrange(grobs = plots, nrow = 2)




##--------------------------------------------------------------------------------------------------------
# merged_cpg_results and matched_tsp_tissue_groups needs to be merged to get cpg info
##--------------------------------------------------------------------------------------------------------

merged_cpg_results_CPG <- lapply(merged_cpg_results, function(df) {
  df = df %>% filter(Category == "High CpG")
  df = df %>% dplyr::select(Header, CpG_Percentage)
  return(df)
})

merged_cpg_results_NONCpg = lapply(merged_cpg_results, function(df) {
  df = df %>% filter(Category == "Low CpG")
  df = df %>% dplyr::select(Header, CpG_Percentage)
  return(df)
}) 

merged_cpg_results_CPG %>% names
# 1] "low_1k_TSS_1k"  "low_2k_TSS_1k"  "low_3k_TSS_1k"  "null_1k_TSS_1k" "null_2k_TSS_1k" "null_3k_TSS_1k" "tenh_1k_TSS_1k" "tenh_2k_TSS_1k" "tenh_3k_TSS_1k" "tsp_1k_TSS_1k"  "tsp_2k_TSS_1k"  "tsp_3k_TSS_1k" "wide_1k_TSS_1k" "wide_2k_TSS_1k" "wide_3k_TSS_1k"

regions <- c("1k_TSS_1k", "2k_TSS_1k", "3k_TSS_1k")
categories <- c("low", "tenh", "wide", "null")

paired_ncpg_combinations <- TSpCpG_otherNonCpG <- list()
paired_cpg_combinations <- TSpNonCpG_otherCpG <- list()

for (region in regions) {
  for (cat in categories) {
    key1 <- paste0("tsp_", region) # maps to merged_cpg_results_CPG: tsp_2k_TSS_1k
    key2 <- paste0(cat, "_", region) # maps to merged_cpg_results_CPG: low_3k_TSS_1k
    name <- paste0(key1, "_vs_", cat)
    
    if (key1 %in% names(merged_cpg_results_NONCpg) && key2 %in% names(merged_cpg_results_NONCpg)) {
      paired_cpg_combinations[[name]] <- rbind(
        merged_cpg_results_CPG[[key1]] %>% mutate(Label = 1),
        merged_cpg_results_CPG[[key2]] %>% mutate(Label = 0))
      
      paired_ncpg_combinations[[name]] <- rbind(
        merged_cpg_results_NONCpg[[key1]] %>% mutate(Label = 1),
        merged_cpg_results_NONCpg[[key2]] %>% mutate(Label = 0))
      
      TSpCpG_otherNonCpG[[name]] <- rbind(
        merged_cpg_results_CPG[[key1]] %>% mutate(Label = 1),
        merged_cpg_results_NONCpg[[key2]] %>% mutate(Label = 0))
      
      TSpNonCpG_otherCpG[[name]] <- rbind(
        merged_cpg_results_NONCpg[[key1]] %>% mutate(Label = 1),
        merged_cpg_results_CPG[[key2]] %>% mutate(Label = 0))
    }}
}

##--------------------------------------------------------------------------------------------------------
# make cpg: different lengthhs like 2k data (1k_TSS_1k)
##--------------------------------------------------------------------------------------------------------

seq_cols <- c("1k_TSS_1k_seq",
              "2k_TSS_1k_seq",
              "3k_TSS_1k_seq")

matched_tsp_tissue_groups1 <- lapply(seq_cols, function(seq_col) {
  # derive start/end column names by replacing "_seq" with "_start"/"_end"
  start_col <- sub("_seq$", "_start", seq_col)
  end_col   <- sub("_seq$", "_end",   seq_col)
  
  # now build the per‐group list
  lapply(matched_tsp_tissue_groups, function(df) {
    df %>%
      mutate(Header = paste0(
        chromosome_name, "_",
        .data[[start_col]],   "_",
        .data[[end_col]]
      )) %>%
      select(Header, all_of(seq_col), transcript, strand, species)
  })
})

# name the outer list by the seq columns
names(matched_tsp_tissue_groups1) <- seq_cols
matched_tsp_tissue_groups1 %>% names
# "1k_TSS_1k_seq"  "2k_TSS_1k_seq"  "3k_TSS_1k_seq" 

############  make the format of below lists like the matched_tsp_tissue_groups1
# paired_ncpg_combinations TSpCpG_otherNonCpG 
# paired_cpg_combinations  TSpNonCpG_otherCpG 

regions <- sub("_seq$", "", seq_cols)
all_pairings <- list(
  TSpCpG_vs_otherNonCpG = TSpCpG_otherNonCpG,
  TSpNonCpG_vs_otherCpG = TSpNonCpG_otherCpG,
  CpG_vs_CpG            = paired_cpg_combinations,      
  NonCpG_vs_NonCpG      = paired_ncpg_combinations  
)

# Build a nested list: one element per region, containing only entries matching that region
nested_pair_lists <- lapply(all_pairings, function(flat_list) {
  # For each region, grab only elements whose names contain that region
  nl <- setNames(
    lapply(regions, function(reg) {
      sel <- grep(reg, names(flat_list), value = TRUE)
      flat_list[sel]
    }),
    paste0(regions, "_seq")
  )
  nl
})




##-----------------------
# final merge
##-----------------------

nested_pair_lists %>% names
# [1] "TSpCpG_vs_otherNonCpG" "TSpNonCpG_vs_otherCpG" "CpG_vs_CpG"            "NonCpG_vs_NonCpG"  

matched_tsp_tissue_groups1 %>% names
nested_pair_lists[[1]] %>% names

seq_keys <- names(matched_tsp_tissue_groups1)
seq_keys

for (nested_list_nm in names(nested_pair_lists)){
  merged_nested <- lapply(seq_keys, function(seq_key, nested_pair_lists_name=nested_list_nm) {
    # sub-lists for this seq
    pair_sub    <- nested_pair_lists[[nested_pair_lists_name]][[seq_key]]
    matched_sub <- matched_tsp_tissue_groups1[[seq_key]]
    
    # extract available categories from each sub-list
    cats_pair    <- sub(".*_vs_", "", names(pair_sub))
    cats_matched <- sub(".*_",   "", names(matched_sub))
    
    # only keep the intersection of categories
    cats <- intersect(cats_pair, cats_matched)
    
    # for each category, merge the corresponding df's
    out <- setNames(lapply(cats, function(cat) {
      # exact names in each sub-list
      pair_name    <- paste0("tsp_", sub("_seq$", "", seq_key), "_vs_", cat)
      match_names  <- grep(paste0("_", cat, "$"), names(matched_sub), value = TRUE)
      
      df1 <- pair_sub[[pair_name]]
      # If multiple matched tables (e.g. multiple tissues for 'null'), bind them
      df2 <- plyr::ldply(matched_sub[match_names])
      colnames(df2)[1] = "experiment"
      
      merge(df1, df2, by = "Header")
    }), cats)
    
    out
  })
  names(merged_nested)=names(matched_tsp_tissue_groups1)
  saveRDS(merged_nested, file.path(save_pth_dna, "rds/cpg_n_cpg/2_processed", paste0(nested_list_nm, ".rds")))
}



##--------------------------------------------------------------------------------------------------------
# 1:1 split then make test train and dev from it
##--------------------------------------------------------------------------------------------------------

# merged_matched_groups
TSpCpG_vs_otherNonCpG = readRDS(file.path(save_pth_dna, "rds/cpg_n_cpg/2_processed", paste0("TSpCpG_vs_otherNonCpG", ".rds")))
TSpNonCpG_vs_otherCpG = readRDS(file.path(save_pth_dna, "rds/cpg_n_cpg/2_processed", paste0("TSpNonCpG_vs_otherCpG", ".rds")))
CpG_vs_CpG = readRDS(file.path(save_pth_dna, "rds/cpg_n_cpg/2_processed", paste0("CpG_vs_CpG", ".rds")))
NonCpG_vs_NonCpG = readRDS(file.path(save_pth_dna, "rds/cpg_n_cpg/2_processed", paste0("NonCpG_vs_NonCpG", ".rds")))

names(TSpCpG_vs_otherNonCpG)=names(TSpNonCpG_vs_otherCpG)=names(CpG_vs_CpG)=names(NonCpG_vs_NonCpG) = c("2000", "3000", "4000")


##--------------------------------------------------------------------------------------------------------
# create split to finetune
##--------------------------------------------------------------------------------------------------------
library(caret)

# now split the low data tsp_brain_low  tsp_liver_low tsp_testis_low based one xperiment column and save all cols and for the test train dev 80:10:10 solit wiill contain Sequence (column ending in _seq) rename it and Label onlyu
out_root="/home/psurana/projects/Davuluri_lab/DNABERT_TSpProm/fine_tune_multispecies_dna2/split_cpg_n_cpg/new"
# dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

set.seed(42)

dirnm="NonCpG_vs_NonCpG"

for (region in names(NonCpG_vs_NonCpG)) {
  for (category in names(NonCpG_vs_NonCpG[[region]])) {
    df_cat <- NonCpG_vs_NonCpG[[region]][[category]]
    seq_col <- grep("_seq$", names(df_cat), value = TRUE)
    
    # iterate over experiments
    for (exp in unique(df_cat$experiment)) {
      df_exp <- df_cat %>% filter(experiment == exp) %>%
                           select(all_of(c(seq_col, "Label"))) %>%
                           rename(Sequence = all_of(seq_col))
      
      counts   <- table(df_sel$Label)
      min_n    <- min(counts)
      df_bal   <- df_sel %>%
        group_by(Label) %>%
        sample_n(min_n) %>%
        ungroup()
      
      # create train/dev/test indices
      train_idx <- createDataPartition(df_bal$Label, p = 0.8, list = FALSE)
      train     <- df_bal[train_idx, ]
      rem       <- df_bal[-train_idx, ]
      
      # now split remainder 50/50 into dev/test = each 10% of original
      dev_idx <- createDataPartition(rem$Label, p = 0.5, list = FALSE)
      dev     <- rem[dev_idx, ]
      test    <- rem[-dev_idx, ]
      
      # make directories
      dir_exp <- file.path(out_root, dirnm, region, exp)
      dir.create(dir_exp, recursive = TRUE, showWarnings = FALSE)
      
      # write CSVs
      write.csv(train, file.path(dir_exp, "train.csv"), row.names = FALSE)
      write.csv(dev,   file.path(dir_exp, "dev.csv"),   row.names = FALSE)
      write.csv(test,  file.path(dir_exp, "test.csv"),  row.names = FALSE)
    }
  }
}

