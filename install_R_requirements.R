# install_R_requirements.R
# Script to install all CRAN and Bioconductor dependencies for TSProm

# ---- 1. Install base package manager ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman", repos = "https://cloud.r-project.org")

# ---- 2. Load and install all packages ----
suppressPackageStartupMessages({
  pacman::p_load(
    optparse, dplyr, stringr, tibble, purrr, data.table, janitor, readr, caret,
    GenomicRanges, IRanges, Biostrings, BSgenome,
    gkmSVM, rtracklayer,
    BSgenome.Hsapiens.UCSC.hg38, 
    BSgenome.Hsapiens.UCSC.hg38.masked,
    BSgenome.Mmusculus.UCSC.mm39
  )
})

# ---- 3. Ensure Bioconductor dependencies ----
BiocManager::install(
  c("GenomicRanges", "IRanges", "Biostrings", "BSgenome",
    "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm39"),
  ask = FALSE, update = FALSE
)
