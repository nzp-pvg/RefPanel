#!/usr/bin/env Rscript
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(glmnet) })
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
mat_file <- arg("--matrix")
meta_file <- arg("--meta")
features_file <- arg("--features")
out <- arg("--out", file.path(data_root, "frozen_model", "C_frozen_model_v1.rds"))
if (any(vapply(list(mat_file, meta_file, features_file), is.null, logical(1)))) {
  stop("Provide archived A inputs with --matrix, --meta and --features. The compact public release distributes the frozen model artifact but not the source matrices.")
}
set.seed(1)
X <- readRDS(mat_file)
meta <- fread(meta_file)
features <- as.character(readRDS(features_file))
if (!all(c("sample_id", "group") %in% names(meta))) stop("A meta requires sample_id/group")
if (is.null(rownames(X)) || is.null(colnames(X))) stop("A matrix needs row/column names")
features <- intersect(features, rownames(X))
if (length(features) < 1000L) stop("Too few frozen features found")
meta <- meta[match(colnames(X), sample_id)]
if (anyNA(meta$group)) stop("A metadata does not align to matrix")
y <- as.integer(meta$group == "HCC")
if (length(unique(y)) != 2L) stop("Training labels need HCC and CTL")
Xf <- t(X[features, , drop = FALSE])
storage.mode(Xf) <- "double"
if (any(!is.finite(Xf))) stop("Non-finite values in training matrix")
cvfit <- cv.glmnet(Xf, y, family = "binomial", alpha = 1, nfolds = 5,
                   type.measure = "auc", grouped = FALSE, standardize = TRUE)
model <- list(model = cvfit, feature_order = features,
              transform = "input matrices must be log2(count_or_signal + 1)",
              glmnet_standardize = TRUE, selected_lambda = cvfit$lambda.1se,
              alpha = 1, nfolds = 5, seed = 1,
              training_samples = colnames(X), training_group = meta$group,
              source_matrix = normalizePath(mat_file), source_features = normalizePath(features_file),
              note = "Frozen on Cohort A only; Cohorts C and D excluded from feature/model fitting")
saveRDS(model, out, compress = "xz")
fwrite(data.table(feature_id = features, feature_order = seq_along(features)), sub("\\.rds$", "_feature_order.tsv", out), sep = "\t")
fwrite(data.table(field = c("n_features", "n_training_samples", "n_HCC", "n_CTL", "lambda_1se", "alpha", "standardize", "transform", "D_used"),
                  value = c(length(features), nrow(meta), sum(y == 1), sum(y == 0), cvfit$lambda.1se, 1, TRUE, "log2(count_or_signal+1)", FALSE)),
       sub("\\.rds$", "_manifest.tsv", out), sep = "\t")
message("[DONE] C frozen model: ", out, " | features=", length(features), " | lambda.1se=", signif(cvfit$lambda.1se, 6))
