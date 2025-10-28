#!/usr/bin/env Rscript

# ======================== #
#  Auto setup & parsing    #
# ======================== #

req_cran <- c("optparse","data.table","dplyr","stringr","tibble","purrr","Biostrings","GenomicRanges","rtracklayer")
req_bioc <- c("gkmSVM","BSgenome","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Hsapiens.UCSC.hg38.masked")
inst <- function(pkgs, bioc=FALSE){
  for(p in pkgs){
    if(!requireNamespace(p, quietly=TRUE)){
      if(bioc){
        if(!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager", repos="https://cloud.r-project.org")
        BiocManager::install(p, ask=FALSE, update=FALSE)
      } else {
        install.packages(p, repos="https://cloud.r-project.org")
      }
    }
    suppressPackageStartupMessages(library(p, character.only=TRUE))
  }
}
inst(req_cran, bioc=FALSE)
inst(req_bioc, bioc=TRUE)

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(stringr); library(tibble)
  library(purrr); library(Biostrings); library(GenomicRanges); library(rtracklayer)
  library(gkmSVM); library(BSgenome.Hsapiens.UCSC.hg38); library(BSgenome.Hsapiens.UCSC.hg38.masked)
})

opt_list <- list(
  make_option("--input_rds", type="character", help="Path to human_seq.rds (or similar)"),
  make_option("--outdir",    type="character", help="Output directory"),
  make_option("--tissues",   type="character", default="spleen,muscle",
              help="Comma-separated tissues (case-insensitive), e.g. 'spleen,muscle'"),
  make_option("--windows",   type="character", default="1k,2k,3k",
              help="Comma-separated window roots that exist in columns like '1k_TSS_1k_start'"),
  make_option("--xfold",     type="integer", default=1, help="genNullSeqs xfold"),
  make_option("--gc_tol",    type="double",  default=0.02, help="GC match tolerance"),
  make_option("--rep_tol",   type="double",  default=0.02, help="Repeat match tolerance"),
  make_option("--len_tol",   type="double",  default=0.02, help="Length match tolerance"),
  make_option("--batch",     type="integer", default=5000, help="genNullSeqs batchsize"),
  make_option("--trials",    type="integer", default=20,   help="genNullSeqs nMaxTrials")
)
opt <- parse_args(OptionParser(option_list=opt_list))
if (is.null(opt$input_rds) || is.null(opt$outdir)) {
  cat("ERROR: --input_rds and --outdir are required.\n"); quit(status=2)
}

input_rds <- normalizePath(opt$input_rds, mustWork=TRUE)
outdir    <- normalizePath(opt$outdir, mustWork=FALSE)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

tissues   <- str_split(opt$tissues, ",")[[1]] %>% trimws() %>% tolower()
wins      <- str_split(opt$windows, ",")[[1]] %>% trimws()

message(">> CONFIG")
message("input_rds : ", input_rds)
message("outdir    : ", outdir)
message("tissues   : ", paste(tissues, collapse=", "))
message("windows   : ", paste(wins, collapse=", "))

# ======================== #
#      Load & prepare      #
# ======================== #

x <- readRDS(input_rds)

# normalize names
colnames(x) <- gsub("\\.", "_", colnames(x))
x <- x %>%
  mutate(
    tissue = tolower(tissue),
    strand = ifelse(strand %in% c("+","-"), strand, ifelse(strand==1,"+","-")),
    chromosome_name = ifelse(startsWith(chromosome_name,"chr"), chromosome_name, paste0("chr", chromosome_name))
  )

# focus on TspTrans
x_tsp <- x %>% filter(tolower(category_final) %in% c("tsp","tsptrans"))
if (nrow(x_tsp)==0) { stop("No TspTrans rows found in input RDS.") }

# build groups list by tissue selection (fallback to tspAll if tissue not present)
present_tissues <- intersect(unique(x_tsp$tissue), tissues)
groups <- list()
if (length(present_tissues) > 0) {
  for (t in present_tissues) groups[[t]] <- x_tsp %>% filter(tissue == t)
} else {
  groups[["tspall"]] <- x_tsp
  message("No requested tissues found; using all TspTrans as 'tspall'.")
}

# ======================== #
# Export BEDs per window   #
# ======================== #

bed_dir <- file.path(outdir, "beds"); dir.create(bed_dir, showWarnings=FALSE, recursive=TRUE)

must_have_cols <- function(df, k) {
  sc <- paste0(k, "_TSS_1k_start")
  ec <- paste0(k, "_TSS_1k_end")
  all(c(sc,ec) %in% colnames(df))
}

write_bed <- function(df, k, name_prefix) {
  sc <- paste0(k, "_TSS_1k_start")
  ec <- paste0(k, "_TSS_1k_end")
  bed_path <- file.path(bed_dir, paste0(name_prefix, "_", k, "TSS1k_pos.bed"))
  bed_df <- df %>%
    select(chromosome_name, !!sym(sc), !!sym(ec), tissue, transcription_start_site, strand) %>%
    as.data.frame()
  write.table(bed_df, bed_path, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  bed_path
}

pos_beds <- list()
for (gname in names(groups)) {
  df <- groups[[gname]]
  for (k in wins) {
    if (!must_have_cols(df, k)) { message("Skip ", gname, " / ", k, " (missing columns)"); next }
    p <- write_bed(df, k, gname)
    pos_beds <- append(pos_beds, p)
    message("Wrote POS BED: ", p)
  }
}

if (length(pos_beds) == 0) stop("No POS BEDs were written (check --windows match your columns).")

# ======================== #
# genNullSeqs (negatives)  #
# ======================== #

# choose genome (currently hg38 masked)
genome_masked <- BSgenome.Hsapiens.UCSC.hg38.masked

for (pos_bed in pos_beds) {
  sample_id <- sub("_pos\\.bed$", "", basename(pos_bed))
  neg_bed   <- file.path(bed_dir, paste0(sample_id, "_neg.bed"))
  pos_fa    <- file.path(bed_dir, paste0(sample_id, "_pos.fa"))
  neg_fa    <- file.path(bed_dir, paste0(sample_id, "_neg.fa"))

  message("genNullSeqs: ", sample_id)
  tryCatch({
    gkmSVM::genNullSeqs(
      inputBedFN        = pos_bed,
      genome            = genome_masked,
      outputBedFN       = neg_bed,
      outputPosFastaFN  = pos_fa,
      outputNegFastaFN  = neg_fa,
      xfold             = opt$xfold,
      repeat_match_tol  = opt$rep_tol,
      GC_match_tol      = opt$gc_tol,
      length_match_tol  = opt$len_tol,
      batchsize         = opt$batch,
      nMaxTrials        = opt$trials
    )
  }, error=function(e){
    message("FAILED genNullSeqs for ", sample_id, ": ", e$message)
  })
}

# ======================== #
# Build labeled outputs    #
# ======================== #

read_fasta_df <- function(fp){
  s <- readDNAStringSet(fp)
  tibble(header = names(s), sequence = as.character(s))
}

# POS (from data frame columns)
make_pos_table <- function(df, k, label_name){
  sc <- paste0(k, "_TSS_1k_start")
  ec <- paste0(k, "_TSS_1k_end")
  qc <- paste0(k, "_TSS_1k_seq")
  if (!all(c(sc,ec,qc) %in% colnames(df))) return(NULL)
  df %>%
    mutate(
      header = paste(chromosome_name, .data[[sc]], .data[[ec]],
                     transcription_start_site, strand, sep = "_"),
      sequence = .data[[qc]],
      sample = paste0(label_name, "_", k, "TSS1k"),
      Label = 1L
    ) %>%
    select(header, sequence, sample, Label)
}

# POS: assemble from groups
pos_all <- list()
for (gname in names(groups)) {
  for (k in wins) {
    tbl <- make_pos_table(groups[[gname]], k, gname)
    if (!is.null(tbl)) pos_all[[length(pos_all)+1]] <- tbl
  }
}
pos_all <- if (length(pos_all)) bind_rows(pos_all) else tibble()

# NEG: read all *_neg.fa under bed_dir
neg_files <- list.files(bed_dir, pattern=".*_neg\\.fa$", full.names=TRUE)
neg_all <- if (length(neg_files)) {
  map_dfr(neg_files, ~{
    d <- read_fasta_df(.x)
    d$sample <- sub("_neg\\.fa$", "", basename(.x))
    d$Label <- 0L
    d
  })
} else tibble()

combined <- bind_rows(pos_all, neg_all)
if (nrow(combined)==0) stop("No combined rows produced.")

# bump header with sample + label for uniqueness
combined <- combined %>%
  mutate(header = paste0(header, "_", sample, "_", Label))

saveRDS(combined, file.path(outdir, "combined_pos_neg.rds"))

# write per-sample FASTA
fasta_dir <- file.path(outdir, "combined_fasta"); dir.create(fasta_dir, showWarnings=FALSE)
write_fa <- function(df, fp){
  dna <- DNAStringSet(df$sequence); names(dna) <- df$header
  writeXStringSet(dna, filepath=fp)
}
combined %>%
  group_by(sample) %>%
  group_walk(~ write_fa(.x, file.path(fasta_dir, paste0(.y$sample, ".fa"))))

message("Done. Outputs:")
message("  RDS  : ", file.path(outdir, "combined_pos_neg.rds"))
message("  FASTA: ", fasta_dir)
