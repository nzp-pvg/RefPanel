#!/usr/bin/env Rscript
rm(list = ls()); gc()
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
features_file <- arg("--features")
manifest_file <- arg("--manifest")
analysis <- arg("--analysis", "primary_cfDNA")
out_prefix <- arg("--out-prefix", "locked_matrix")
if (is.null(features_file) || is.null(manifest_file)) stop("Required: --features FILE --manifest FILE")
features <- fread(features_file)
if (!all(c("chr", "start", "end", "feature_id") %in% names(features))) stop("Feature TSV lacks chr/start/end/feature_id")
manifest <- fread(manifest_file)
if (!analysis %in% names(manifest)) stop("Manifest lacks analysis flag: ", analysis)
use <- manifest[get(analysis) == TRUE]
if (!nrow(use)) stop("No included samples")
key <- paste(features$chr, features$start, features$end, sep = "\t")

read_locked <- function(path) {
  x <- fread(path, header = FALSE, select = 1:4, col.names = c("chr", "start", "end", "value"))
  idx <- match(key, paste(x$chr, x$start, x$end, sep = "\t"))
  if (anyNA(idx)) stop("Locked coordinates missing from ", path, ": ", sum(is.na(idx)))
  as.numeric(x$value[idx])
}
mat <- vapply(use$path, read_locked, numeric(nrow(features)))
rownames(mat) <- features$feature_id
colnames(mat) <- use$sample_id
saveRDS(mat, paste0(out_prefix, ".rds"), compress = "xz")
fwrite(data.table(feature_id = rownames(mat), mat), paste0(out_prefix, ".tsv.gz"), sep = "\t")
fwrite(use, paste0(out_prefix, "_samples.tsv"), sep = "\t")
qc <- data.table(n_features = nrow(mat), n_samples = ncol(mat), any_NA = anyNA(mat),
                 min_value = min(mat), max_value = max(mat), analysis = analysis)
fwrite(qc, paste0(out_prefix, "_qc.tsv"), sep = "\t")
message("[DONE] ", nrow(mat), " locked features x ", ncol(mat), " samples")
