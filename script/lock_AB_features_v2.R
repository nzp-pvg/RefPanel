#!/usr/bin/env Rscript
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
a_bed <- arg("--a-stable-bed")
b_meta <- arg("--b-meta")
out_dir <- arg("--out", "out_AB_lock_v2")
ic_quantile <- as.numeric(arg("--ic-quantile", "0.20"))
pseudocount <- as.numeric(arg("--pseudocount", "0.1"))
if (is.null(a_bed) || is.null(b_meta)) stop("Required: --a-stable-bed FILE --b-meta FILE")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bed <- fread(a_bed, header = FALSE)
if (ncol(bed) < 4L) stop("A stable BED requires four columns")
setnames(bed, 1:4, c("chr", "start", "end", "feature_id"))
if (any(bed$end - bed$start + 1L != 300L)) stop("Primary lock is fixed at 300 bp; non-300 feature found")
meta <- fread(b_meta)
if (!all(c("key", "type", "file") %in% names(meta))) stop("B meta requires key, type, file")
meta[, type_std := fifelse(grepl("IC", type), "IC", fifelse(grepl("(^|_)M($|_)", type), "M", NA_character_))]
if (anyNA(meta$type_std)) stop("Unrecognized B type")
if (any(meta[, .N, by = .(key, type_std)]$N != 1L)) stop("B keys must have exactly one M and one IC")
meta_dir <- dirname(normalizePath(b_meta))
meta[, file := ifelse(grepl("^(~|/|[A-Za-z]:)", file), file, file.path(meta_dir, file))]
if (any(!file.exists(meta$file))) stop("Missing B avg files; paths in B meta must be executable")

read_mean0 <- function(path, wanted) {
  x <- fread(path, header = FALSE, select = c(1L, 5L))
  setnames(x, c("feature_id", "mean0"))
  x[match(wanted, feature_id), mean0]
}
ic_rows <- meta[type_std == "IC"]
IC <- vapply(ic_rows$file, read_mean0, numeric(nrow(bed)), wanted = bed$feature_id)
colnames(IC) <- ic_rows$key
if (anyNA(IC)) stop("At least one A feature is absent from a B IC file")
ic_median <- apply(IC, 1L, median)
threshold <- unname(quantile(ic_median, ic_quantile, type = 8))
keep <- ic_median > threshold
locked <- bed[keep]
locked[, `:=`(scale_bp = 300L, A_definition = "CTL_signaleligible_lowest10pct_CV",
              B_IC_median_mean0 = ic_median[keep], B_IC_gate_quantile = ic_quantile,
              B_IC_threshold_mean0 = threshold, ratio_pseudocount = pseudocount)]
fwrite(locked, file.path(out_dir, "AB_locked_features_v2.tsv"), sep = "\t")
fwrite(locked[, .(chr, start, end, feature_id)], file.path(out_dir, "AB_locked_features_v2.bed"),
       sep = "\t", col.names = FALSE)
manifest <- data.table(field = c("feature_selection_cohorts", "target_cohorts_excluded", "scale_bp",
                                 "A_stability_rule", "B_IC_gate_quantile", "B_IC_threshold_mean0",
                                 "ratio_pseudocount", "n_A_candidates", "n_locked"),
                       value = c("A,B", "C,D", 300, "signal eligible then lowest 10% CTL CV",
                                 ic_quantile, threshold, pseudocount, nrow(bed), nrow(locked)))
fwrite(manifest, file.path(out_dir, "AB_locked_features_manifest_v2.tsv"), sep = "\t")
message("[DONE] Locked ", nrow(locked), " A/B-only 300-bp features")
