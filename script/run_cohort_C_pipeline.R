# Cohort C pipeline: preprocess + stability analysis (public version)
# This script is designed to run from the project root using relative paths.

rm(list = ls())
gc()

# ---- 0) configuration ----
# Project root is assumed to be the current working directory.
# If you run this from elsewhere, set PROJECT_ROOT manually.
PROJECT_ROOT <- getwd()
DATA_DIR <- file.path(PROJECT_ROOT, "data", "cohort_C")
A_BED <- file.path(PROJECT_ROOT, "data", "cohort_B", "A_bins_300bp_hg19_blacklist_filtered.bed")

# Output files (written under DATA_DIR)
OUT_BINS <- file.path(DATA_DIR, "C_bins_intersection.txt")
OUT_META <- file.path(DATA_DIR, "C_meta.tsv")
OUT_CF   <- file.path(DATA_DIR, "C_cfDNA_log2_matrix.rds")
OUT_PBL  <- file.path(DATA_DIR, "C_pbl_log2_matrix.rds")
OUT_NCC  <- file.path(DATA_DIR, "C_NCC_log2_matrix.rds")

# ---- 1) dependencies ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(matrixStats)
})

# ---- 2) helpers ----
logi <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
}

# Load an RData file and return the object named 'matri'.
load_matrix_rdata <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  objs <- ls(e)
  stopifnot(length(objs) >= 1)
  if (!("matri" %in% objs)) stop("No object named 'matri' in: ", path)
  get("matri", envir = e)
}

# Normalize bin IDs: chr1.750301.750600 -> chr1_750301_750600
normalize_bins <- function(x) {
  rn <- rownames(x)
  rownames(x) <- gsub("\\.", "_", rn)
  x
}

# Get upper-triangle values of a correlation matrix (no diagonal).
upper_vals <- function(C) {
  C[upper.tri(C, diag = FALSE)]
}

# Summarize Spearman correlations using bin subsampling.
spearman_summary_subsample <- function(mat_log2, n_row = 50000, B = 3, seed = 1) {
  set.seed(seed)
  p <- ncol(mat_log2)
  stopifnot(p >= 2)

  summarise_cor <- function(C) {
    v <- upper_vals(C)
    data.frame(
      n_row = nrow(C),
      n_samples = p,
      median_spearman = median(v, na.rm = TRUE),
      IQR_spearman = IQR(v, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  out <- vector("list", B)
  nr_total <- nrow(mat_log2)

  for (b in seq_len(B)) {
    idx <- sample.int(nr_total, min(n_row, nr_total))
    sub <- mat_log2[idx, , drop = FALSE]
    C <- suppressWarnings(cor(sub, method = "spearman", use = "pairwise.complete.obs"))

    s <- summarise_cor(C)
    s$n_row <- nrow(sub)
    s$replicate <- b
    out[[b]] <- s
  }

  do.call(rbind, out)
}

# ---- 3) inputs ----
if (!dir.exists(DATA_DIR)) stop("Cannot find data directory: ", DATA_DIR)
if (!file.exists(A_BED)) stop("Cannot find A bins bed: ", A_BED)

setwd(DATA_DIR)

logi("[READ] A bins bed: ", A_BED)
bed <- read.delim(A_BED, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
if (ncol(bed) < 4) stop("BED must have >=4 columns: chr,start,end,name")
A_bins <- bed[[4]]
rm(bed)
gc()
logi("[OK] A bins loaded: ", length(A_bins))

# ---- 4) load matrices (raw counts) ----
logi("[READ] cfDNA counts RData...")
cf <- load_matrix_rdata("cfDNA_medip_counts_matrix_blacklist_filtered.RData")
logi("[OK] cfDNA dim: ", paste(dim(cf), collapse = " x "))

logi("[READ] pbl counts RData...")
pbl <- load_matrix_rdata("pbl_medip_counts_matrix_blacklist_filtered.RData")
logi("[OK] pbl dim: ", paste(dim(pbl), collapse = " x "))

logi("[READ] NCC counts RData...")
ncc <- load_matrix_rdata("NCC_cfdna_medip_matrix.RData")
logi("[OK] NCC dim: ", paste(dim(ncc), collapse = " x "))

# ---- 5) harmonize rownames ----
cf  <- normalize_bins(cf)
pbl <- normalize_bins(pbl)
ncc <- normalize_bins(ncc)

# ---- 6) intersect bins with A bins ----
common_bins <- Reduce(intersect, list(rownames(cf), rownames(pbl), rownames(ncc), A_bins))
logi("[OK] common bins (C matrices ∩ A bins): ", length(common_bins))

cf  <- cf[common_bins, , drop = FALSE]
pbl <- pbl[common_bins, , drop = FALSE]
ncc <- ncc[common_bins, , drop = FALSE]

# ---- 7) log2(count+1) transform ----
logi("[DO] log2(count+1) transform...")
cf_log2  <- log2(cf + 1)
pbl_log2 <- log2(pbl + 1)
ncc_log2 <- log2(ncc + 1)

rm(cf, pbl, ncc)
gc()

# ---- 8) metadata ----
logi("[READ] Clinical Data.csv")
clin <- read.csv("Clinical Data.csv", stringsAsFactors = FALSE)

meta_cf <- tibble(sample_id = colnames(cf_log2)) %>%
  mutate(group = "SCLC") %>%
  left_join(clin %>% rename(sample_id = ID), by = "sample_id")

meta_ncc <- tibble(sample_id = colnames(ncc_log2)) %>%
  mutate(group = "NCC")

meta <- bind_rows(meta_cf, meta_ncc)

# ---- 9) save outputs ----
writeLines(common_bins, OUT_BINS)
write.table(meta, OUT_META, sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(cf_log2,  OUT_CF,  compress = FALSE)
saveRDS(pbl_log2, OUT_PBL, compress = FALSE)
saveRDS(ncc_log2, OUT_NCC, compress = FALSE)

logi("[DONE] Saved:")
logi("  - ", OUT_CF)
logi("  - ", OUT_PBL)
logi("  - ", OUT_NCC)
logi("  - ", OUT_META)
logi("  - ", OUT_BINS)
logi("[INFO] final dim: ", nrow(cf_log2), " bins | cfDNA=", ncol(cf_log2),
     " | pbl=", ncol(pbl_log2), " | NCC=", ncol(ncc_log2))

# ---- 10) NCC stability analysis (optional) ----
# This block quantifies sample-to-sample Spearman stability with bin subsampling.

ncc_res_50k  <- spearman_summary_subsample(ncc_log2, n_row = 50000,  B = 3, seed = 1)
ncc_res_200k <- spearman_summary_subsample(ncc_log2, n_row = 200000, B = 3, seed = 2)
ncc_res_500k <- spearman_summary_subsample(ncc_log2, n_row = 500000, B = 3, seed = 3)

ncc_res <- rbind(ncc_res_50k, ncc_res_200k, ncc_res_500k)
write.table(ncc_res, file.path(DATA_DIR, "C_NCC_spearman_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ncc_sum <- ncc_res %>%
  group_by(n_row, n_samples) %>%
  summarise(
    median_of_median = median(median_spearman),
    IQR_of_median    = IQR(median_spearman),
    mean_median      = mean(median_spearman),
    sd_median        = sd(median_spearman),
    median_IQR       = median(IQR_spearman),
    .groups = "drop"
  )

write.table(ncc_sum, file.path(DATA_DIR, "C_NCC_spearman_sensitivity_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 11) cfDNA stability analysis (optional) ----
cf_res_50k  <- spearman_summary_subsample(cf_log2,  n_row = 50000,  B = 2, seed = 11)
cf_res_200k <- spearman_summary_subsample(cf_log2,  n_row = 200000, B = 2, seed = 12)
cf_res_500k <- spearman_summary_subsample(cf_log2,  n_row = 500000, B = 2, seed = 13)

cf_res <- rbind(cf_res_50k, cf_res_200k, cf_res_500k)
write.table(cf_res, file.path(DATA_DIR, "C_cfDNA_spearman_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
