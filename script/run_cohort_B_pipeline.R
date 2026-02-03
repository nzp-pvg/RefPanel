#!/usr/bin/env Rscript
rm(list = ls()); gc()

### making meta data

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(tidyr)
})

# Working directory: cohort_B
# Expect bigWigAverageOverBed outputs under avg_dir (each file is an avg.txt)
avg_dir <- "avg"
stopifnot(dir.exists(avg_dir))

files <- list.files(avg_dir, full.names = TRUE)
files <- files[grepl("\\.txt$|\\.tsv$|\\.avg$", files, ignore.case = TRUE)]
stopifnot(length(files) > 0)

# Parse visit (A/C) and type (M/IC) from filename.
# This is intentionally permissive: it searches anywhere in the basename.
parse_one <- function(f) {
  bn <- basename(f)

  visit <- NA_character_
  if (str_detect(bn, "(^|[^A-Za-z0-9])A([^A-Za-z0-9]|$)")) visit <- "A"
  if (str_detect(bn, "(^|[^A-Za-z0-9])C([^A-Za-z0-9]|$)")) visit <- "C"

  type <- NA_character_
  if (str_detect(bn, "(^|[^A-Za-z0-9])IC([^A-Za-z0-9]|$)")) type <- "IC"
  if (str_detect(bn, "(^|[^A-Za-z0-9])M([^A-Za-z0-9]|$)"))  type <- "M"

  # A conservative "subject id" extraction:
  # remove trailing parts starting from visit/type tokens; keep the left trunk as subject id
  subject <- bn
  subject <- str_replace(subject, "(^|_)(A|C)(_|\\.).*$", "")  # if visit exists
  subject <- str_replace(subject, "(^|_)(IC|M)(_|\\.).*$", "") # if type exists
  subject <- str_replace(subject, "\\.txt$|\\.tsv$|\\.avg$", "")

  tibble(file = f, fname = bn, subject = subject, visit = visit, type = type)
}

meta <- bind_rows(lapply(files, parse_one))

# Keep only visit A/C and type M/IC rows
meta_ac <- meta %>%
  filter(visit %in% c("A","C")) %>%
  filter(type %in% c("M","IC")) %>%
  mutate(
    key = paste0(subject, "_", visit)
  ) %>%
  select(key, subject, visit, type, file)

# Report parsing failures explicitly
bad <- meta %>% filter(!(visit %in% c("A","C")) | !(type %in% c("M","IC")))
if (nrow(bad) > 0) {
  cat("\n[WARN] Some files could not be parsed into (visit,type). Examples:\n")
  print(head(bad %>% select(fname, visit, type), 30))
  cat("\n[WARN] Those files are ignored. Fix filename tokens if needed.\n\n")
}

stopifnot(nrow(meta_ac) > 0)

# Save the raw meta used by build_matrix
write_tsv(meta_ac, "B_meta_AC.tsv")

# ---- Clean version: standardize type and enforce paired M+IC for each key ----
meta_clean <- read_tsv("B_meta_AC.tsv", show_col_types = FALSE) %>%
  mutate(
    type_std = case_when(
      str_detect(type, "IC") ~ "IC",
      str_detect(type, "M")  ~ "M",
      TRUE ~ NA_character_
    )
  )

stopifnot(all(!is.na(meta_clean$type_std)))

chk <- meta_clean %>%
  count(key, type_std, name = "n") %>%
  pivot_wider(names_from = type_std, values_from = n, values_fill = 0)

bad_pair <- chk %>% filter(M != 1 | IC != 1)
if (nrow(bad_pair) > 0) {
  cat("\n[ERROR] Some keys are not perfectly paired (need exactly 1 M and 1 IC):\n")
  print(bad_pair)
  stop("Fix duplicates/missing files, then regenerate meta.")
}

write_tsv(meta_clean, "B_meta_AC_clean.tsv")

cat("\n[DONE] Wrote:\n")
cat("  - B_meta_AC.tsv\n")
cat("  - B_meta_AC_clean.tsv\n")
cat("Keys:", n_distinct(meta_clean$key), " | rows:", nrow(meta_clean), "\n")



suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
})

# ----------------------------
# Config
# ----------------------------
data_dir <- "data"
out_dir  <- "out"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_in  <- file.path(data_dir, "B_meta_AC.tsv")
stopifnot(file.exists(meta_in))

visits_keep <- c("A", "C")

eps_ratio <- 1e-3      # for ratio log2((M+eps)/(IC+eps))
pc_ratio  <- 0.1       # for sensitivity ratio log2((M+pc)/(IC+pc))
topN      <- 5000

q_list <- seq(0.05, 0.50, by = 0.05)
main_q <- 0.20

set.seed(1)

# ----------------------------
# Helpers
# ----------------------------
ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logi <- function(...) message(sprintf("[%s] %s", ts(), paste0(..., collapse = "")))

as_relpath <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^(~|/|[A-Za-z]:)", x), x, file.path(data_dir, x))
}

read_mean0 <- function(f) {
  # bigWigAverageOverBed output: name size covered sum mean0 mean
  dt <- data.table::fread(f, header = FALSE)
  stopifnot(ncol(dt) >= 5)
  v <- dt[[5]]
  names(v) <- dt[[1]]
  v
}

topN_overlap <- function(x, y, topN = 5000) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 10) return(NA_real_)
  topN2 <- min(topN, length(x))
  ix <- order(x, decreasing = TRUE)[seq_len(topN2)]
  iy <- order(y, decreasing = TRUE)[seq_len(topN2)]
  length(intersect(ix, iy)) / topN2
}

summ_median_iqr <- function(x) {
  x <- x[is.finite(x)]
  c(median = median(x, na.rm = TRUE), IQR = IQR(x, na.rm = TRUE), n = length(x))
}

# ----------------------------
# B01: meta clean + paired check + build matrices
# ----------------------------
logi("B01 start: read meta")
meta <- read.delim(meta_in, sep = "\t", stringsAsFactors = FALSE) |> tibble::as_tibble()

stopifnot(all(c("key", "type", "visit", "file") %in% colnames(meta)))

meta <- meta |>
  dplyr::mutate(
    type_std = dplyr::case_when(
      grepl("(^IC\\b|\\b_IC\\b)", type) ~ "IC",
      grepl("(^M\\b|\\b_M\\b)",  type) ~ "M",
      TRUE ~ NA_character_
    ),
    visit = as.character(visit),
    file  = as_relpath(file)
  )

stopifnot(all(!is.na(meta$type_std)))
stopifnot(all(meta$visit %in% visits_keep))
stopifnot(all(file.exists(meta$file)))

meta <- meta |> dplyr::filter(visit %in% visits_keep)

chk <- meta |>
  dplyr::count(key, type_std, name = "n") |>
  tidyr::pivot_wider(names_from = type_std, values_from = n, values_fill = 0)

bad <- chk |> dplyr::filter(M != 1 | IC != 1)
if (nrow(bad) > 0) {
  print(bad)
  stop("Some keys are missing exactly one M and one IC.")
}
logi("B01 paired check OK")

meta_out <- file.path(out_dir, "B_meta_AC_clean.tsv")
write.table(meta, meta_out, sep = "\t", quote = FALSE, row.names = FALSE)
logi("B01 saved: ", meta_out)

logi("B01 build tile universe from first file")
x0 <- read_mean0(meta$file[1])
tiles <- names(x0)
stopifnot(length(tiles) > 10)

keys <- sort(unique(meta$key))
M_mat  <- matrix(NA_real_, nrow = length(tiles), ncol = length(keys),
                 dimnames = list(tiles, keys))
IC_mat <- matrix(NA_real_, nrow = length(tiles), ncol = length(keys),
                 dimnames = list(tiles, keys))

for (i in seq_len(nrow(meta))) {
  f <- meta$file[i]
  k <- meta$key[i]
  t <- meta$type_std[i]
  logi("READ ", k, " ", t, " | ", basename(f))

  v <- read_mean0(f)
  v <- v[tiles]

  if (t == "M")  M_mat[, k]  <- v
  if (t == "IC") IC_mat[, k] <- v
}

stopifnot(identical(dim(M_mat), dim(IC_mat)))
stopifnot(identical(rownames(M_mat), rownames(IC_mat)))
stopifnot(identical(colnames(M_mat), colnames(IC_mat)))

stopifnot(all(is.finite(colSums(M_mat, na.rm = TRUE))))
stopifnot(all(is.finite(colSums(IC_mat, na.rm = TRUE))))

saveRDS(M_mat,  file.path(out_dir, "B_M_mean0_matrix_AC.rds"))
saveRDS(IC_mat, file.path(out_dir, "B_IC_mean0_matrix_AC.rds"))
logi("B01 saved: B_M_mean0_matrix_AC.rds, B_IC_mean0_matrix_AC.rds")

# ----------------------------
# B015: ratio main
# ----------------------------
ratio_log2 <- log2((M_mat + eps_ratio) / (IC_mat + eps_ratio))
saveRDS(ratio_log2, file.path(out_dir, "B_ratio_log2_matrix_AC.rds"))
logi("B015 saved: B_ratio_log2_matrix_AC.rds")

# also save M log2 for convenience
M_log2 <- log2(M_mat + 1)
saveRDS(M_log2, file.path(out_dir, "B_M_log2_matrix_AC.rds"))
logi("B015 saved: B_M_log2_matrix_AC.rds")

# ----------------------------
# B02: consistency MAIN
# ----------------------------
logi("B02 start: consistency metrics")
keys <- colnames(M_log2)
stopifnot(identical(keys, colnames(ratio_log2)))

spearman_bykey <- vapply(keys, function(k) {
  x <- M_log2[, k]
  y <- ratio_log2[, k]
  suppressWarnings(stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}, numeric(1))

top_overlap_bykey <- vapply(keys, function(k) {
  x <- M_log2[, k]
  y <- ratio_log2[, k]
  topN_overlap(x, y, topN = topN)
}, numeric(1))

res_cons <- data.frame(
  key = keys,
  spearman = spearman_bykey,
  topN = topN,
  topN_overlap_frac = top_overlap_bykey,
  stringsAsFactors = FALSE
)

f_cons <- file.path(out_dir, "B_consistency_bykey.tsv")
write.table(res_cons, f_cons, sep = "\t", quote = FALSE, row.names = FALSE)

summ_cons <- data.frame(
  metric = c("spearman(log2(M+1), log2(M/IC))", sprintf("top%d overlap fraction", topN)),
  median = c(median(res_cons$spearman, na.rm = TRUE),
             median(res_cons$topN_overlap_frac, na.rm = TRUE)),
  IQR    = c(IQR(res_cons$spearman, na.rm = TRUE),
             IQR(res_cons$topN_overlap_frac, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

f_summ <- file.path(out_dir, "B_consistency_summary.tsv")
write.table(summ_cons, f_summ, sep = "\t", quote = FALSE, row.names = FALSE)

p1 <- ggplot(res_cons, aes(x = spearman)) +
  geom_histogram(bins = 20) +
  theme_bw() +
  labs(x = "Spearman per key", y = "Count",
       title = "Cohort B: input only vs ratio consistency")

ggsave(file.path(out_dir, "B_fig_consistency_spearman.png"), p1, width = 7, height = 4, dpi = 150)

p2 <- ggplot(res_cons, aes(x = topN_overlap_frac)) +
  geom_histogram(bins = 20) +
  theme_bw() +
  labs(x = sprintf("Top %d overlap fraction per key", topN), y = "Count",
       title = "Cohort B: top feature concordance")

ggsave(file.path(out_dir, "B_fig_consistency_topN_overlap.png"), p2, width = 7, height = 4, dpi = 150)

logi("B02 saved: B_consistency_bykey.tsv, B_consistency_summary.tsv, 2 figures")

# ----------------------------
# B03: sensitivity and scatter
# ----------------------------
logi("B03 start: sensitivity sweep")

stopifnot(identical(dim(M_mat), dim(IC_mat)))
ic_tile_med <- matrixStats::rowMedians(IC_mat, na.rm = TRUE)

all_res <- list()

for (qq in q_list) {
  thr <- as.numeric(stats::quantile(ic_tile_med, probs = qq, na.rm = TRUE))
  keep <- ic_tile_med > thr
  keep_n <- sum(keep)

  logi("B03 q=", sprintf("%.2f", qq),
       " thr=", signif(thr, 5),
       " keep=", keep_n, "/", length(keep))

  M_log <- log2(M_mat[keep, , drop = FALSE] + 1)
  ratio <- log2((M_mat[keep, , drop = FALSE] + pc_ratio) / (IC_mat[keep, , drop = FALSE] + pc_ratio))

  spearman_vec <- vapply(keys, function(k) {
    suppressWarnings(stats::cor(M_log[, k], ratio[, k], method = "spearman", use = "pairwise.complete.obs"))
  }, numeric(1))

  overlap_vec <- vapply(keys, function(k) {
    topN_overlap(M_log[, k], ratio[, k], topN = topN)
  }, numeric(1))

  df <- data.frame(
    key = keys,
    ic_filter_quantile = qq,
    pc = pc_ratio,
    topN = topN,
    kept_tiles = keep_n,
    spearman = spearman_vec,
    topN_overlap = overlap_vec,
    stringsAsFactors = FALSE
  )

  out_q <- file.path(out_dir, sprintf("B_sensitivity_bykey_q%02d.tsv", round(qq * 100)))
  write.table(df, out_q, sep = "\t", quote = FALSE, row.names = FALSE)

  all_res[[as.character(qq)]] <- df
  rm(M_log, ratio); gc()
}

res_sens <- dplyr::bind_rows(all_res)
saveRDS(res_sens, file.path(out_dir, "B_sensitivity_res.rds"), compress = "xz")
write.table(res_sens, file.path(out_dir, "B_sensitivity_res.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

sum_tbl <- res_sens |>
  dplyr::group_by(ic_filter_quantile) |>
  dplyr::summarise(
    spearman_median = median(spearman, na.rm = TRUE),
    spearman_IQR    = IQR(spearman, na.rm = TRUE),
    overlap_median  = median(topN_overlap, na.rm = TRUE),
    overlap_IQR     = IQR(topN_overlap, na.rm = TRUE),
    n_keys          = dplyr::n(),
    .groups = "drop"
  )

write.table(sum_tbl, file.path(out_dir, "B_sensitivity_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

df_main <- res_sens |> dplyr::filter(abs(ic_filter_quantile - main_q) < 1e-12)

p_scatter <- ggplot(df_main, aes(x = spearman, y = topN_overlap)) +
  geom_point(size = 2, alpha = 0.85) +
  labs(
    title = "Cohort B: concordance per key",
    subtitle = sprintf("IC filter quantile=%.2f, pc=%.1f, TopN=%d", main_q, pc_ratio, topN),
    x = "Spearman: log2(M+1) vs log2((M+pc)/(IC+pc))",
    y = sprintf("Top %d overlap fraction", topN)
  ) +
  theme_bw()

ggsave(file.path(out_dir, "B_fig_scatter_spearman_vs_overlap_MAIN.png"),
       p_scatter, width = 10, height = 6, dpi = 150)

sum_long <- sum_tbl |>
  dplyr::select(ic_filter_quantile, spearman_median, spearman_IQR, overlap_median, overlap_IQR) |>
  tidyr::pivot_longer(cols = c(spearman_median, overlap_median),
                      names_to = "metric", values_to = "median") |>
  dplyr::mutate(
    IQR = dplyr::if_else(
      metric == "spearman_median",
      sum_tbl$spearman_IQR[match(ic_filter_quantile, sum_tbl$ic_filter_quantile)],
      sum_tbl$overlap_IQR[match(ic_filter_quantile, sum_tbl$ic_filter_quantile)]
    ),
    metric = dplyr::if_else(
      metric == "spearman_median",
      "Spearman (median, IQR)",
      sprintf("Top %d overlap (median, IQR)", topN)
    )
  )

p_sens <- ggplot(sum_long, aes(x = ic_filter_quantile, y = median)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = median - IQR / 2, ymax = median + IQR / 2), width = 0.01) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  labs(
    title = "Cohort B: sensitivity to IC tile gating",
    subtitle = sprintf("pc=%.1f; error bars show IQR/2", pc_ratio),
    x = "IC filter quantile (remove lowest tiles by median(IC))",
    y = "Median metric value"
  ) +
  theme_bw()

ggsave(file.path(out_dir, "B_fig_sensitivity_ICquantile.png"),
       p_sens, width = 10, height = 8, dpi = 150)

logi("B03 saved: per q TSV, combined TSV RDS, summary TSV, 2 figures")

# ----------------------------
# Console markdown summary
# ----------------------------
cat("\n## Cohort B: consistency (MAIN)\n\n")
cat("| Metric | Median | IQR |\n|---|---:|---:|\n")
cat(sprintf("| %s | %.4f | %.4f |\n", summ_cons$metric[1], summ_cons$median[1], summ_cons$IQR[1]))
cat(sprintf("| %s | %.6f | %.6f |\n", summ_cons$metric[2], summ_cons$median[2], summ_cons$IQR[2]))

cat("\n## Cohort B: sensitivity summary\n\n")
cat("| IC_filter_quantile | spearman_median | spearman_IQR | overlap_median | overlap_IQR | n_keys |\n")
cat("|---:|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(sum_tbl))) {
  cat(sprintf("| %.2f | %.4f | %.4f | %.4f | %.4f | %d |\n",
              sum_tbl$ic_filter_quantile[i],
              sum_tbl$spearman_median[i],
              sum_tbl$spearman_IQR[i],
              sum_tbl$overlap_median[i],
              sum_tbl$overlap_IQR[i],
              sum_tbl$n_keys[i]))
}

cat("\n[DONE] outputs are in: ", out_dir, "\n", sep = "")
