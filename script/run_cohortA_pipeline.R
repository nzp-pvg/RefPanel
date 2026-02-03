## Raw data can be downloaded from 🔽🔽🔽
https://zenodo.org/records/11251606


#!/usr/bin/env Rscript
rm(list=ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(matrixStats)
})

# ----------------------------
# Config
# ----------------------------
data_dir <- "data"
out_dir  <- "out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts_rdata <- file.path(data_dir, "vCount_n236.rds")
sample_file  <- file.path(data_dir, "sample.rds")
black_rdata  <- file.path(data_dir, "black_bin_v2.RData")

scales <- c(300, 600, 900, 1200)   # bp
top_frac <- 0.1
B_boot <- 25
boot_frac <- 0.8
set.seed(1)

# ----------------------------
# Helpers
# ----------------------------
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logi <- function(...) message(sprintf("[%s] %s", ts(), paste0(..., collapse="")))

jaccard <- function(a, b) {
  a <- unique(a); b <- unique(b)
  length(intersect(a, b)) / length(union(a, b))
}

pairwise_spearman_median <- function(mat) {
  C <- suppressWarnings(cor(mat, method = "spearman"))
  C[lower.tri(C, diag = TRUE)] <- NA_real_
  median(C, na.rm = TRUE)
}

stable_tiles_by_cv <- function(mat, top_frac = 0.1) {
  mu <- rowMeans(mat)
  sdv <- matrixStats::rowSds(mat)
  cv <- sdv / pmax(mu, 1e-6)
  thr <- as.numeric(stats::quantile(cv, probs = top_frac, na.rm = TRUE))
  which(cv <= thr)
}

bootstrap_jaccard <- function(mat, B = 25, frac = 0.8, top_frac = 0.1) {
  n <- ncol(mat)
  k <- max(2, floor(frac * n))
  sets <- vector("list", B)

  for (b in seq_len(B)) {
    cols <- sample.int(n, k)
    idx <- stable_tiles_by_cv(mat[, cols, drop = FALSE], top_frac = top_frac)
    sets[[b]] <- idx
  }

  jac <- c()
  for (i in 1:(B - 1)) for (j in (i + 1):B) {
    a <- sets[[i]]; b <- sets[[j]]
    jac <- c(jac, length(intersect(a, b)) / length(union(a, b)))
  }
  mean(jac)
}

parse_bins_from_rownames <- function(rn) {
  df <- data.frame(bin_id = rn, stringsAsFactors = FALSE)
  tmp <- stringr::str_split_fixed(df$bin_id, "_", 3)
  df$chr   <- tmp[, 1]
  df$start <- as.integer(tmp[, 2])
  df$end   <- as.integer(tmp[, 3])
  df
}

write_bed4 <- function(df_bins, out_file) {
  bed4 <- df_bins[, c("chr","start","end","bin_id")]
  write.table(bed4, out_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}

merge_bins <- function(counts300, bins300, k = 2) {
  o <- order(bins300$chr, bins300$start)
  bins300 <- bins300[o, ]
  counts300 <- counts300[o, , drop = FALSE]

  grp <- ave(seq_len(nrow(bins300)), bins300$chr, FUN=function(ix) ((ix-1) %/% k) + 1)

  out_bins <- bins300 %>%
    dplyr::mutate(grp = grp) %>%
    dplyr::group_by(chr, grp) %>%
    dplyr::summarise(start = min(start), end = max(end), .groups="drop") %>%
    dplyr::mutate(bin_id = paste0(chr, "_", start, "_", end))

  out_counts <- rowsum(counts300, group = paste0(bins300$chr, "___", grp), reorder = FALSE)
  rownames(out_counts) <- out_bins$bin_id

  list(counts = out_counts, bins = out_bins)
}

# ----------------------------
# 1) Load inputs
# ----------------------------
stopifnot(file.exists(counts_rdata), file.exists(sample_file), file.exists(black_rdata))
logi("Loading counts: ", counts_rdata)

e <- new.env(parent = emptyenv())
load(counts_rdata, envir = e)
obj_names <- ls(e)
stopifnot(length(obj_names) == 1)
counts_all <- as.matrix(e[[obj_names]])
rm(e); gc()

logi("Loading sample sheet: ", sample_file)
samples <- data.table::fread(sample_file, data.table = FALSE)
stopifnot(all(c("Sample_ID","Group") %in% colnames(samples)))

logi("Loading blacklist: ", black_rdata)
e2 <- new.env(parent = emptyenv())
load(black_rdata, envir = e2)
stopifnot("black_bin" %in% ls(e2))
black_bin <- e2[["black_bin"]]
rm(e2); gc()
stopifnot("black_bin" %in% colnames(black_bin))

stopifnot(!is.null(rownames(counts_all)), !is.null(colnames(counts_all)))

logi("counts dim = ", paste(dim(counts_all), collapse=" x "))
logi("samples dim = ", paste(dim(samples), collapse=" x "))
logi("blacklist n = ", nrow(black_bin))

# ----------------------------
# 2) Subset CTL and apply blacklist
# ----------------------------
sid <- as.character(samples$Sample_ID)
cid <- colnames(counts_all)
common <- intersect(sid, cid)
logi("common samples = ", length(common))

samples_use <- samples %>% dplyr::filter(Sample_ID %in% common)

ctl_ids <- samples_use %>%
  dplyr::filter(Group == "CTL") %>%
  dplyr::pull(Sample_ID) %>%
  as.character()

ctl_ids <- setdiff(ctl_ids, "N19")
stopifnot(length(ctl_ids) >= 3)
logi("CTL n = ", length(ctl_ids))

counts_ctl <- counts_all[, ctl_ids, drop = FALSE]

black_set <- unique(as.character(black_bin$black_bin))
keep <- !(rownames(counts_ctl) %in% black_set)
counts_ctl_blk_300 <- counts_ctl[keep, , drop = FALSE]
logi("counts (300bp, CTL, blacklist) dim = ", paste(dim(counts_ctl_blk_300), collapse=" x "))

# save meta
write.table(
  samples_use %>% dplyr::filter(Sample_ID %in% ctl_ids),
  file = file.path(out_dir, "A_meta_CTL.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ----------------------------
# 3) Build 300bp bed and multi-scale bins
# ----------------------------
bins300 <- parse_bins_from_rownames(rownames(counts_ctl_blk_300))
stopifnot(all(is.finite(bins300$start)), all(is.finite(bins300$end)))
stopifnot(all(bins300$end - bins300$start == 299))

# 300bp outputs
saveRDS(counts_ctl_blk_300, file.path(out_dir, "A_CTL_counts_raw_300bp_blacklist_filtered.rds"))
write_bed4(bins300, file.path(out_dir, "A_bins_300bp_blacklist_filtered.bed"))

# merged scales
for (bp in scales[scales != 300]) {
  k <- bp / 300
  res <- merge_bins(counts_ctl_blk_300, bins300, k = k)
  saveRDS(res$counts, file.path(out_dir, sprintf("A_CTL_counts_raw_%sbp_blacklist_filtered.rds", bp)))
  write_bed4(res$bins, file.path(out_dir, sprintf("A_bins_%sbp_blacklist_filtered.bed", bp)))
  logi("saved ", bp, "bp counts+bed | dim=", paste(dim(res$counts), collapse=" x "))
}

# ----------------------------
# 4) Metrics + stable tiles per scale
# ----------------------------
metrics_rows <- list()
stable_sets <- list()

for (bp in scales) {
  f_counts <- file.path(out_dir, sprintf("A_CTL_counts_raw_%sbp_blacklist_filtered.rds", bp))
  stopifnot(file.exists(f_counts))

  logi("metrics start: ", bp, "bp")
  m <- as.matrix(readRDS(f_counts))

  lcpm <- log2(m + 1)

  med_cor <- pairwise_spearman_median(lcpm)

  mu <- rowMeans(m)
  cv <- matrixStats::rowSds(m) / pmax(mu, 1e-6)
  madv <- matrixStats::rowMads(lcpm)

  stable_idx <- stable_tiles_by_cv(m, top_frac = top_frac)
  stable_prop <- length(stable_idx) / nrow(m)

  jac <- bootstrap_jaccard(m, B = B_boot, frac = boot_frac, top_frac = top_frac)

  stable_bins <- rownames(m)[stable_idx]
  out_stable_rds <- file.path(out_dir, sprintf("A_stable_bins_top10pctCV_%sbp.rds", bp))
  saveRDS(stable_bins, out_stable_rds)

  metrics_rows[[as.character(bp)]] <- data.frame(
    bin_size_bp = bp,
    n_tiles = nrow(m),
    n_samples = ncol(m),
    median_pairwise_spearman_log = med_cor,
    median_CV_raw = median(cv, na.rm = TRUE),
    median_MAD_log2 = median(madv, na.rm = TRUE),
    stable_tiles_top10pct_prop = stable_prop,
    bootstrap_stabletiles_jaccard_mean = jac,
    stringsAsFactors = FALSE
  )

  stable_sets[[as.character(bp)]] <- stable_bins

  rm(m, lcpm, mu, cv, madv, stable_idx, stable_bins); gc()
  logi("metrics done: ", bp, "bp")
}

metrics_tbl <- dplyr::bind_rows(metrics_rows)
write.table(metrics_tbl, file.path(out_dir, "A_metrics_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# markdown print
md <- metrics_tbl
md[] <- lapply(md, function(x) if (is.numeric(x)) signif(x, 4) else x)
cat("\n## Cohort A stability metrics (CTL, blacklist-filtered)\n\n")
cat("| bin_size_bp | n_tiles | n_samples | median_pairwise_spearman_log | median_CV_raw | median_MAD_log2 | stable_tiles_top10pct_prop | bootstrap_stabletiles_jaccard_mean |\n")
cat("|---:|---:|---:|---:|---:|---:|---:|---:|\n")
for (r in seq_len(nrow(md))) cat("|", paste(md[r, ], collapse=" | "), "|\n")

# ----------------------------
# 5) Cross-scale Jaccard matrix (stable bins)
# ----------------------------
J <- matrix(NA_real_, nrow = length(scales), ncol = length(scales),
            dimnames = list(as.character(scales), as.character(scales)))

for (i in seq_along(scales)) {
  for (j in seq_along(scales)) {
    J[i, j] <- jaccard(stable_sets[[as.character(scales[i])]],
                       stable_sets[[as.character(scales[j])]])
  }
}

write.table(
  cbind(bin_size_bp = rownames(J), as.data.frame(J, check.names = FALSE)),
  file = file.path(out_dir, "A_stablebins_jaccard_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n[DONE] wrote out/\n",
    "A_CTL_counts_raw_*bp_blacklist_filtered.rds\n",
    "A_bins_*bp_blacklist_filtered.bed\n",
    "A_stable_bins_top10pctCV_*bp.rds\n",
    "A_metrics_summary.tsv\n",
    "A_stablebins_jaccard_matrix.tsv\n", sep = "")

# ----------------------------
# 6) Export stable bins as BED per scale
# ----------------------------
for (bp in scales) {
  bed_file <- file.path(out_dir, sprintf("A_bins_%sbp_blacklist_filtered.bed", bp))
  stable_rds <- file.path(out_dir, sprintf("A_stable_bins_top10pctCV_%sbp.rds", bp))
  out_bed <- file.path(out_dir, sprintf("A_stable_bins_top10pctCV_%sbp.bed", bp))

  stopifnot(file.exists(bed_file), file.exists(stable_rds))

  bed <- data.table::fread(bed_file, header = FALSE)
  stopifnot(ncol(bed) >= 4)
  colnames(bed)[1:4] <- c("chr","start","end","name")

  stable_ids <- readRDS(stable_rds)
  bed2 <- bed[name %in% stable_ids, .(chr, start, end, name)]
  data.table::fwrite(bed2, out_bed, sep = "\t", col.names = FALSE)

  logi("stable bed written: ", out_bed, " | n=", nrow(bed2))
}
