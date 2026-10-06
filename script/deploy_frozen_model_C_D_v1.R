#!/usr/bin/env Rscript
rm(list = ls()); gc()
suppressPackageStartupMessages({ library(data.table); library(glmnet) })
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
args <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default = NULL) { i <- match(flag, args); if (is.na(i)) default else args[[i + 1L]] }
model_file <- arg("--model", file.path(data_root, "frozen_model", "C_frozen_model_v1.rds"))
out <- arg("--out", file.path(data_root, "deployment"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
model <- readRDS(model_file); fit <- model$model; features <- model$feature_order
auc_rank <- function(prob, y) { n1 <- sum(y == 1); n0 <- sum(y == 0); if (!n1 || !n0) return(NA_real_); r <- rank(prob); (sum(r[y == 1]) - n1*(n1+1)/2)/(n1*n0) }
score <- function(mat, ids, groups, cohort, material) {
  miss <- setdiff(features, rownames(mat)); if (length(miss)) stop(cohort, " missing frozen features: ", length(miss))
  x <- t(mat[features, , drop = FALSE]); if (any(!is.finite(x))) stop(cohort, " non-finite values")
  prob <- as.numeric(predict(fit, newx = x, s = model$selected_lambda, type = "response"))
  y <- as.integer(groups == "HCC" | groups == "HNSCC" | groups == "SCLC")
  data.table(cohort = cohort, material = material, sample_id = ids, group = groups, frozen_score = prob,
             positive_label = y, auc = auc_rank(prob, y))
}

# Cohort C: optional because the public compact release contains final scores, not its full matrix.
c_res <- NULL
c_cf <- arg("--c-cfdna")
c_ncc <- arg("--c-ncc")
c_meta <- arg("--c-meta")
if (!is.null(c_cf) || !is.null(c_ncc) || !is.null(c_meta)) {
  if (any(vapply(list(c_cf, c_ncc, c_meta), is.null, logical(1)))) stop("Provide --c-cfdna, --c-ncc and --c-meta together")
  c <- cbind(readRDS(c_cf), readRDS(c_ncc)); cm <- fread(c_meta)[match(colnames(c), sample_id)]
  if (anyNA(cm$group)) stop("C metadata mismatch")
  c_res <- score(c, colnames(c), as.character(cm$group), "C", "cfDNA")
  fwrite(c_res, file.path(out, "C_frozen_scores_v1.tsv"), sep = "\t")
}

# Cohort D: locked compact cfDNA and PBL matrices.
ddir <- file.path(data_root, "cohort_D", "derived_v1")
d_manifest <- rbind(fread(file.path(ddir, "D_cfDNA_model20k_v1_samples.tsv")), fread(file.path(ddir, "D_PBL_model20k_v1_samples.tsv")), fill=TRUE)
d_res <- rbind(score(readRDS(file.path(ddir, "D_cfDNA_model20k_v1.rds")), d_manifest[primary_cfDNA == TRUE, sample_id],
                     ifelse(d_manifest[primary_cfDNA == TRUE, status] == "patient", "HNSCC", "healthy"), "D", "cfDNA"),
              score(readRDS(file.path(ddir, "D_PBL_model20k_v1.rds")), d_manifest[secondary_PBL == TRUE, sample_id],
                    ifelse(d_manifest[secondary_PBL == TRUE, status] == "patient", "HNSCC", "healthy"), "D", "PBL"))
fwrite(d_res, file.path(out, "D_frozen_scores_v1.tsv"), sep = "\t")
sum_c <- if (is.null(c_res)) NULL else c_res[, .(cohort = unique(cohort), material = unique(material), n = .N, n_positive = sum(positive_label), n_negative = sum(!positive_label), auc = unique(auc))]
sum_d <- d_res[, .(cohort = unique(cohort), n = .N, n_positive = sum(positive_label), n_negative = sum(!positive_label), auc = unique(auc)), by = material]
fwrite(rbindlist(Filter(Negate(is.null), list(sum_c, sum_d)), fill = TRUE),
       file.path(out, "frozen_deployment_summary_v1.tsv"), sep = "\t")
message("[DONE] frozen scores written: ", normalizePath(out))
