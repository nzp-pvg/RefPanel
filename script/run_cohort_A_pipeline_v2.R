#!/usr/bin/env Rscript
rm(list = ls()); gc()

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) return(default)
  if (hit == length(args)) stop("Missing value after ", flag)
  args[[hit + 1L]]
}

counts_file <- arg("--counts", "data/vCount_n236.rds")
sample_file <- arg("--samples", "data/sample.rds")
black_file <- arg("--blacklist", "data/black_bin_v2.RData")
out_dir <- arg("--out", "out_A_v2")
top_frac <- as.numeric(arg("--stable-frac", "0.10"))
min_mean_count <- as.numeric(arg("--min-mean-count", "1"))
n_subsamples <- as.integer(arg("--n-subsamples", "25"))
sample_frac <- as.numeric(arg("--sample-frac", "0.80"))
seed <- as.integer(arg("--seed", "1"))
scales <- c(300L, 600L, 900L, 1200L)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_single_object <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    x <- try(readRDS(path), silent = TRUE)
    if (!inherits(x, "try-error")) return(x)
  }
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  nm <- ls(e)
  if (length(nm) == 1L) return(e[[nm]])
  # Common legacy A RData contains both the count matrix and black_bin.
  if ("black_bin" %in% nm && any(grepl("count|vCount", nm, ignore.case = TRUE))) {
    return(e[[nm[grepl("count|vCount", nm, ignore.case = TRUE)][1L]]])
  }
  stop("Expected one object in ", path, "; found: ", paste(nm, collapse = ", "))
}

read_samples <- function(path) {
  x <- try(read_single_object(path), silent = TRUE)
  if (!inherits(x, "try-error") && is.data.frame(x)) return(x)
  x <- try(fread(path, data.table = FALSE), silent = TRUE)
  if (!inherits(x, "try-error") && is.data.frame(x)) return(x)
  stop("Sample input must be a readable RDS/RData or delimited table: ", path)
}

parse_bins <- function(ids) {
  z <- tstrsplit(ids, "_", fixed = TRUE)
  if (length(z) != 3L) stop("Bin IDs must be chr_start_end")
  data.table(chr = z[[1]], start = as.integer(z[[2]]), end = as.integer(z[[3]]), bin_id = ids)
}

make_scale <- function(counts, bins, bp, black_ids) {
  k <- bp %/% 300L
  setorder(bins, chr, start, end)
  counts <- counts[match(bins$bin_id, rownames(counts)), , drop = FALSE]
  bins[, grid_index := start %/% 300L]
  bins[, group_index := grid_index %/% k]
  bins[, group_key := paste(chr, group_index, sep = "___")]
  qc <- bins[, .(
    n_parts = .N,
    out_start = min(start),
    out_end = max(end),
    consecutive = all(diff(start) == 300L),
    component_blacklisted = any(bin_id %in% black_ids)
  ), by = group_key]
  qc[, valid := n_parts == k & consecutive & (out_end - out_start + 1L) == bp & !component_blacklisted]
  valid_keys <- qc[valid == TRUE, group_key]
  out <- rowsum(counts, bins$group_key, reorder = FALSE)
  out <- out[valid_keys, , drop = FALSE]
  valid_qc <- qc[match(valid_keys, group_key)]
  valid_qc[, chr := sub("___.*$", "", valid_keys)]
  ids <- paste(valid_qc$chr, valid_qc$out_start, valid_qc$out_end, sep = "_")
  rownames(out) <- ids
  bed <- data.table(chr = valid_qc$chr, start = valid_qc$out_start, end = valid_qc$out_end, bin_id = ids)
  list(counts = out, bed = bed, qc = qc)
}

stable_ids <- function(m, frac, min_mean) {
  mu <- rowMeans(m)
  cv <- rowSds(m) / pmax(mu, .Machine$double.eps)
  eligible <- is.finite(cv) & mu >= min_mean
  if (sum(eligible) < 10L) stop("Too few signal-eligible bins; check --min-mean-count")
  threshold <- unname(quantile(cv[eligible], frac, na.rm = TRUE, type = 8))
  list(ids = rownames(m)[eligible & cv <= threshold], mean = mu, cv = cv,
       eligible = eligible, cv_threshold = threshold)
}

subsample_jaccard <- function(m, B, sample_frac, frac, min_mean, seed) {
  set.seed(seed)
  k <- max(2L, floor(ncol(m) * sample_frac))
  sets <- replicate(B, stable_ids(m[, sample.int(ncol(m), k), drop = FALSE], frac, min_mean)$ids,
                    simplify = FALSE)
  vals <- unlist(lapply(seq_len(B - 1L), function(i) vapply((i + 1L):B, function(j) {
    length(intersect(sets[[i]], sets[[j]])) / length(union(sets[[i]], sets[[j]]))
  }, numeric(1))))
  mean(vals)
}

atomic_300_ids <- function(bed) {
  unlist(Map(function(chr, start, end) {
    starts <- seq(start, end - 299L, by = 300L)
    paste(chr, starts, starts + 299L, sep = "_")
  },
             bed$chr, bed$start, bed$end), use.names = FALSE)
}

counts <- as.matrix(read_single_object(counts_file))
samples <- read_samples(sample_file)
be <- new.env(parent = emptyenv()); load(black_file, envir = be)
if (!"black_bin" %in% ls(be)) stop("Blacklist RData lacks object black_bin: ", black_file)
black <- be[["black_bin"]]
if (!all(c("Sample_ID", "Group") %in% names(samples))) stop("Samples require Sample_ID and Group columns")
if (!"black_bin" %in% names(black)) stop("Blacklist requires black_bin column")
if (is.null(rownames(counts)) || is.null(colnames(counts))) stop("Counts require row and column names")
bins300 <- parse_bins(rownames(counts))
if (anyNA(bins300$start) || anyNA(bins300$end) || any(bins300$end - bins300$start + 1L != 300L))
  stop("Counts are not a valid 300-bp chr_start_end grid")
ctl <- setdiff(as.character(samples$Sample_ID[samples$Group == "CTL"]), "N19")
ctl <- ctl[ctl %in% colnames(counts)]
if (length(ctl) < 3L) stop("Fewer than three matched CTL samples")
counts <- counts[, ctl, drop = FALSE]
black_ids <- unique(as.character(black$black_bin))

metric <- list(); stable_beds <- list(); atomic_sets <- list()
for (bp in scales) {
  obj <- make_scale(counts, copy(bins300), bp, black_ids)
  s <- stable_ids(obj$counts, top_frac, min_mean_count)
  bed_stable <- obj$bed[match(s$ids, bin_id)]
  fwrite(obj$bed, file.path(out_dir, sprintf("A_bins_%dbp_contiguous_blacklist_filtered.bed", bp)),
         sep = "\t", col.names = FALSE)
  fwrite(bed_stable, file.path(out_dir, sprintf("A_stable_signaleligible_top10pctCV_%dbp.bed", bp)),
         sep = "\t", col.names = FALSE)
  saveRDS(s$ids, file.path(out_dir, sprintf("A_stable_signaleligible_top10pctCV_%dbp.rds", bp)))
  metric[[as.character(bp)]] <- data.table(
    bin_size_bp = bp, n_bins = nrow(obj$counts), n_CTL = ncol(obj$counts),
    min_mean_count = min_mean_count, n_signal_eligible = sum(s$eligible),
    stable_fraction_among_eligible = length(s$ids) / sum(s$eligible), cv_threshold = s$cv_threshold,
    subsample_stable_set_jaccard_mean = subsample_jaccard(obj$counts, n_subsamples, sample_frac,
                                                          top_frac, min_mean_count, seed))
  stable_beds[[as.character(bp)]] <- bed_stable
  atomic_sets[[as.character(bp)]] <- unique(atomic_300_ids(bed_stable))
}
fwrite(rbindlist(metric), file.path(out_dir, "A_metrics_v2.tsv"), sep = "\t")
J <- outer(as.character(scales), as.character(scales), Vectorize(function(a, b) {
  length(intersect(atomic_sets[[a]], atomic_sets[[b]])) / length(union(atomic_sets[[a]], atomic_sets[[b]]))
}))
dimnames(J) <- list(scales, scales)
fwrite(data.table(bin_size_bp = rownames(J), as.data.frame(J, check.names = FALSE)),
       file.path(out_dir, "A_cross_scale_atomic300_jaccard.tsv"), sep = "\t")
manifest <- data.table(parameter = c("counts_file", "sample_file", "blacklist_file", "scales_bp",
                                     "stable_fraction", "min_mean_count", "n_subsamples", "sample_fraction", "seed"),
                       value = c(normalizePath(counts_file), normalizePath(sample_file), normalizePath(black_file),
                                 paste(scales, collapse = ","), top_frac, min_mean_count, n_subsamples, sample_frac, seed))
fwrite(manifest, file.path(out_dir, "A_run_manifest_v2.tsv"), sep = "\t")
message("[DONE] Cohort A v2 outputs: ", normalizePath(out_dir))
