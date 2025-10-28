# =========================== #
# TSprom modular pipeline     #
# =========================== #

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(c("GenomicRanges", "IRanges", "Biostrings", 
                       "BSgenome", "BSgenome.Hsapiens.UCSC.hg38", 
                       "BSgenome.Mmusculus.UCSC.mm39"), 
                     ask = FALSE, update = FALSE)

suppressPackageStartupMessages({
  pacman::p_load(
    data.table, dplyr, stringr, janitor, readr, caret,
    GenomicRanges, IRanges, Biostrings, BSgenome,
    BSgenome.Hsapiens.UCSC.hg38, BSgenome.Mmusculus.UCSC.mm39
  )
})

# ---------- CONFIG ---------- #

get_env <- function(var, default = NULL) {
  val <- Sys.getenv(var)
  if (val == "" && !is.null(default)) return(default)
  return(val)
}

config <- list(
  save_base = get_env("SAVE_DIR", "/tmp/tsprom_run"),
  human_csv = get_env("HUMAN_CSV", ""),
  mouse_rds = get_env("MOUSE_RDS", ""),
  nonprom_rds = get_env("NONPROM_RDS", ""),
  species = strsplit(get_env("SPECIES", "human"), ",")[[1]],
  tsp_site_tissues = strsplit(get_env("TISSUES", "brain,testis,liver"), ",")[[1]],
  train_p = as.numeric(get_env("TRAIN_P", 0.7)),
  dev_p = as.numeric(get_env("DEV_P", 0.15)),
  tss_min_separation = as.integer(get_env("TSS_MIN_SEP", 0)),
  drop_with_N = TRUE,
  windows = list(
    "1k_TSS_1k"  = c(-1000, 1000),
    "2k_TSS_1k"  = c(-2000, 1000),
    "3k_TSS_1k"  = c(-3000, 1000)
  ),
  tsp_label = "TspTrans",
  other_labels = c("NullTrans","LowTrans","TenhTrans","WideTrans"),
  seed = 20251028
)

message("✅ Config loaded from environment:")
print(config)

# ---------- HELPERS ---------- #

ensure_dir <- function(...) { dir.create(file.path(...), recursive = TRUE, showWarnings = FALSE) }
nzchar_ <- function(x) !is.na(x) & x != ""

normalize_species_df <- function(df, is_human = TRUE, tsp_label = "TspTrans") {
  # Required cols (names may already match your files)
  expect <- c("chromosome_name","transcript_start","transcript_end",
              "transcription_start_site","strand","tissue","category_final")
  missing <- setdiff(expect, names(df))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse=", "))

  # Keep valid chromosomes
  if (is_human) {
    df <- df %>% filter(chromosome_name %in% c(as.character(1:22),"X","Y"))
  } else {
    df <- df %>% filter(chromosome_name %in% c(as.character(1:19),"X","Y"))
  }

  df %>%
    mutate(
      species = ifelse(is_human, "human","mouse"),
      chromosome_name = paste0("chr", chromosome_name),
      strand = ifelse(strand == 1, "+", "-"),
      category_final = dplyr::case_when(
        category_final == tsp_label ~ "tsp",
        category_final == "NullTrans" ~ "null",
        category_final == "LowTrans"  ~ "low",
        category_final == "TenhTrans" ~ "tenh",
        category_final == "WideTrans" ~ "wide",
        TRUE ~ tolower(category_final)
      ),
      tissue = tolower(tissue)
    ) %>%
    distinct(transcription_start_site, chromosome_name, strand, .keep_all = TRUE)
}

build_gr_for_window <- function(df, up_dn) {
  up <- up_dn[1]; dn <- up_dn[2]
  starts <- ifelse(df$strand == "+",
                   df$transcription_start_site + up,
                   df$transcription_start_site + up) # genomic BED coords don't flip with strand
  ends   <- ifelse(df$strand == "+",
                   df$transcription_start_site + dn,
                   df$transcription_start_site + dn)
  GRanges(seqnames = df$chromosome_name,
          ranges = IRanges(start = pmax(1, as.integer(starts)),
                           end   = pmax(1, as.integer(ends))),
          strand = df$strand)
}

extract_window_sequences <- function(df, windows, genome) {
  out <- df
  for (nm in names(windows)) {
    gr <- build_gr_for_window(df, windows[[nm]])
    seqs <- getSeq(genome, gr)
    out[[paste0(nm, "_seq")]] <- as.character(seqs)
  }
  out
}

drop_and_balance <- function(df, drop_with_N = TRUE) {
  if (drop_with_N) {
    seq_cols <- grep("_seq$", names(df), value = TRUE)
    for (sc in seq_cols) df <- df[!grepl("N", df[[sc]], fixed = TRUE), ]
  }
  # balance by Label
  if (!"Label" %in% names(df)) return(df)
  set.seed(config$seed)
  n0 <- sum(df$Label == 0); n1 <- sum(df$Label == 1)
  if (n0 == 0 || n1 == 0) return(df)
  n   <- min(n0, n1)
  rbind(
    df %>% filter(Label == 1) %>% sample_n(n),
    df %>% filter(Label == 0) %>% sample_n(n)
  )
}

split_70_15_15 <- function(df, train_p = .7, dev_p = .15, seed = 1) {
  stopifnot(all(c("Sequence","Label") %in% names(df)))
  set.seed(seed)
  idx_tr <- createDataPartition(df$Label, p = train_p, list = FALSE)
  train  <- df[idx_tr, ]
  rest   <- df[-idx_tr, ]
  dev_p2 <- dev_p / (1 - train_p)
  idx_de <- createDataPartition(rest$Label, p = dev_p2, list = FALSE)
  dev    <- rest[idx_de, ]
  test   <- rest[-idx_de, ]
  list(train=train, dev=dev, test=test)
}

write_split_csvs <- function(splits, outdir) {
  ensure_dir(outdir)
  write.csv(splits$train[, c("Sequence","Label")], file.path(outdir,"train.csv"), row.names = FALSE)
  write.csv(splits$dev[,   c("Sequence","Label")], file.path(outdir,"dev.csv"),   row.names = FALSE)
  write.csv(splits$test[,  c("Sequence","Label")], file.path(outdir,"test.csv"),  row.names = FALSE)
}

make_header <- function(df) {
  paste(df$chromosome_name, df$transcription_start_site, df$strand,
        ifelse(nzchar_("transcript") && "transcript" %in% names(df), df$transcript, ""),
        df$tissue, df$species, sep = "|")
}

# ---------- MAIN ---------- #

run_pipeline <- function(config) {
  message(">> Starting pipeline…")
  set.seed(config$seed)

  save_rds_dir   <- file.path(config$save_base, "rds");           ensure_dir(save_rds_dir)
  fasta_dir_base <- file.path(config$save_base, "fasta");         ensure_dir(fasta_dir_base)
  csv_dir_base   <- file.path(config$save_base, "datasets");      ensure_dir(csv_dir_base)

  # --- Load & normalize ---
  if ("human" %in% config$species) {
    dh <- fread(config$human_csv, data.table = FALSE)
    # rename transcript id if needed
    if (!"transcript" %in% names(dh) && "ensembl_transcript_id" %in% names(dh)) {
      dh$transcript <- dh$ensembl_transcript_id
    }
    human <- normalize_species_df(dh, is_human=TRUE, tsp_label=config$tsp_label)
  } else human <- NULL

  if ("mouse" %in% config$species) {
    dm <- readRDS(config$mouse_rds)
    if (!"transcript" %in% names(dm) && "ensembl_transcript_id" %in% names(dm)) {
      dm$transcript <- dm$ensembl_transcript_id
    }
    mouse <- normalize_species_df(dm, is_human=FALSE, tsp_label=config$tsp_label)
  } else mouse <- NULL

  df_all <- bind_rows(human, mouse)
  if (!nrow(df_all)) stop("No rows after loading inputs.")

  # Optional min spacing between TSS (by chr)
  if (config$tss_min_separation > 0) {
    df_all <- df_all %>%
      arrange(chromosome_name, transcription_start_site) %>%
      group_by(chromosome_name) %>%
      mutate(TSS_diff = c(NA, diff(transcription_start_site))) %>%
      filter(is.na(TSS_diff) | TSS_diff >= config$tss_min_separation) %>%
      ungroup()
  }

  # --- Extract sequences for each species separately (faster/safer) ---
  seqd <- list()
  if (!is.null(human) && nrow(human)) {
    message("Extracting human windows…")
    seqd$human <- extract_window_sequences(human, config$windows, BSgenome.Hsapiens.UCSC.hg38)
  }
  if (!is.null(mouse) && nrow(mouse)) {
    message("Extracting mouse windows…")
    seqd$mouse <- extract_window_sequences(mouse, config$windows, BSgenome.Mmusculus.UCSC.mm39)
  }
  df_seq <- bind_rows(seqd)

  saveRDS(df_seq, file.path(save_rds_dir, "seq_extracted_all.rds"))

  # --- Define groups & labels ---
  # Positive = tsp; Negative = others
  df_pos <- df_seq %>% filter(category_final == "tsp") %>% mutate(Label = 1L)
  df_neg <- df_seq %>% filter(category_final %in% tolower(config$other_labels)) %>% mutate(Label = 0L)

  # 1) TSP (all tissues/species) vs REST categories
  tspAll_rest <- bind_rows(df_pos, df_neg)

  # 2) TSP site (per selected tissue) vs REST categories
  site_groups <- lapply(config$tsp_site_tissues, function(tt) {
    tt <- tolower(tt)
    pos_site <- df_pos %>% filter(tissue == tt)
    if (!nrow(pos_site)) return(NULL)
    bind_rows(pos_site, df_neg)
  })
  names(site_groups) <- paste0("tsp_", tolower(config$tsp_site_tissues), "_vs_rest")
  site_groups <- Filter(Negate(is.null), site_groups)

  # 3) TSP (human|mouse|both) vs Non-Prom (if provided)
  nonprom_groups <- list()
  if (!is.null(config$nonprom_rds) && file.exists(config$nonprom_rds)) {
    bg <- readRDS(config$nonprom_rds) # expects columns named like the *_seq windows
    # harmonize columns: carry only seq columns and make Label=0
    seq_cols <- grep("_seq$", names(df_seq), value = TRUE)
    for (nm in seq_cols) if (!nm %in% names(bg)) {
      warning("Non-promoter RDS missing column: ", nm, " (skipping bg pairing for this window)")
    }
    bg$Label <- 0L
    # Build pos from df_pos (all tsp, all tissues/species)
    nonprom_groups[["tspAll_vs_nonProm"]] <- list(pos = df_pos, bg = bg)
    # Per-tissue TSP vs non-prom
    for (tt in config$tsp_site_tissues) {
      tt <- tolower(tt)
      pos_site <- df_pos %>% filter(tissue == tt)
      if (nrow(pos_site)) {
        nonprom_groups[[paste0("tsp", tt, "_vs_nonProm")]] <- list(pos = pos_site, bg = bg)
      }
    }
  }

  # --- Writer helpers for datasets/FASTA ---
  write_group_outputs <- function(df_group, group_name) {
    if (is.null(df_group) || !nrow(df_group)) return(invisible(NULL))
    df_group <- df_group %>% select(chromosome_name, transcription_start_site, strand, transcript, tissue, species,
                                    category_final, Label, dplyr::ends_with("_seq"))

    # Drop & balance
    df_group <- drop_and_balance(df_group, config$drop_with_N)

    # Loop windows
    for (win in names(config$windows)) {
      seq_col <- paste0(win, "_seq")
      if (!seq_col %in% names(df_group)) next
      sub <- df_group %>% filter(nzchar_(!!as.name(seq_col)))
      if (!nrow(sub)) next

      # Build {Sequence, Label}
      ds <- sub %>%
        transmute(
          Sequence = .data[[seq_col]],
          Label = Label
        )

      # Split
      splits <- split_70_15_15(ds, train_p = config$train_p, dev_p = config$dev_p, seed = config$seed)

      # Save
      outdir <- file.path(csv_dir_base, group_name, win)
      write_split_csvs(splits, outdir)

      # Optional FASTA (headers are informative)
      headers <- make_header(sub)
      dna <- DNAStringSet(sub[[seq_col]])
      names(dna) <- headers
      fasta_path <- file.path(fasta_dir_base, group_name)
      ensure_dir(fasta_path)
      writeXStringSet(dna, filepath = file.path(fasta_path, paste0(group_name, "_", win, ".fasta")))
    }
  }

  # --- Emit (1) tspAll_vs_rest ---
  write_group_outputs(tspAll_rest, "tspAll_vs_rest")

  # --- Emit (2) per-tissue tsp vs rest ---
  for (nm in names(site_groups)) {
    write_group_outputs(site_groups[[nm]], nm)
  }

  # --- Emit (3) tsp vs non-prom (when available) ---
  if (length(nonprom_groups)) {
    for (nm in names(nonprom_groups)) {
      pair <- nonprom_groups[[nm]]
      # Align columns: keep only windows present in both
      seq_cols_pos <- grep("_seq$", names(pair$pos), value = TRUE)
      seq_cols_bg  <- grep("_seq$", names(pair$bg), value = TRUE)
      shared       <- intersect(seq_cols_pos, seq_cols_bg)
      if (!length(shared)) next

      # Construct a synthetic df with shared windows and labels
      dfp <- pair$pos %>% mutate(Label = 1L) %>% select(Label, all_of(shared))
      dfb <- pair$bg  %>% mutate(Label = 0L) %>% select(Label, all_of(shared))

      # We’ll attach minimal meta for FASTA headers (synthetic where missing)
      # Use pos meta where possible, fallback for bg
      meta_pos <- pair$pos %>% select(chromosome_name, transcription_start_site, strand, transcript, tissue, species) %>%
        mutate(across(everything(), ~replace_na(.x, "")))
      if (!nrow(meta_pos)) next
      # For simplicity, recycle the first meta row for bg headers (or leave blanks)
      meta_bg <- meta_pos[rep(1, nrow(dfb)), , drop = FALSE]

      combo <- list(pos = cbind(meta_pos, dfp), bg = cbind(meta_bg, dfb))
      df_combo <- bind_rows(combo$pos, combo$bg)

      # Balance / drop Ns
      df_combo <- drop_and_balance(df_combo, config$drop_with_N)

      for (win in shared) {
        sub <- df_combo %>% filter(nzchar_(!!as.name(win)))
        if (!nrow(sub)) next
        ds <- sub %>% transmute(Sequence = .data[[win]], Label = Label)
        splits <- split_70_15_15(ds, train_p = config$train_p, dev_p = config$dev_p, seed = config$seed)
        outdir <- file.path(csv_dir_base, nm, gsub("_seq$", "", win))
        write_split_csvs(splits, outdir)

        # FASTA
        headers <- make_header(sub)
        dna <- DNAStringSet(sub[[win]])
        names(dna) <- headers
        ensure_dir(fasta_dir_base, nm)
        writeXStringSet(dna, filepath = file.path(fasta_dir_base, nm, paste0(nm, "_", gsub("_seq$","",win), ".fasta")))
      }
    }
  }

  message(">> Done. Outputs:")
  message("   RDS:   ", save_rds_dir)
  message("   CSVs:  ", csv_dir_base)
  message("   FASTA: ", file.path(fasta_dir_base))
}

# -------- RUN -------- #
run_pipeline(config)
