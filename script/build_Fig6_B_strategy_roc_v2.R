#!/usr/bin/env Rscript

# Cohort D ROC visualization from frozen inputs. No D-specific feature selection,
# scaling, fitting, or threshold optimization is performed here.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."), mustWork = TRUE)
data_root <- file.path(root, "intermediate_data")
derived <- file.path(data_root, "cohort_D", "derived_v1")
deploy <- file.path(data_root, "deployment")
fig_dir <- file.path(root, "figures", "Fig_6")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cf_path <- file.path(derived, "D_cfDNA_model20k_v1.rds")
pbl_path <- file.path(derived, "D_PBL_model20k_v1.rds")
locked_path <- file.path(deploy, "D_frozen_scores_v1.tsv")
bench_path <- file.path(deploy, "D_strategy_benchmark_v1.tsv")
out_scores <- file.path(deploy, "D_strategy_scores_reconstructed_v2.tsv")
required <- c(cf_path, pbl_path, locked_path, bench_path)
if (!all(file.exists(required))) stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))

auc_rank <- function(score, y) {
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1L) / 2) / (n1 * n0)
}
roc_table <- function(d) {
  ord <- order(d$score, decreasing = TRUE)
  y <- d$positive_label[ord]
  data.frame(material = d$material[1], strategy = d$strategy[1], auc = d$auc[1],
             fpr = c(0, cumsum(y == 0L) / sum(y == 0L)),
             tpr = c(0, cumsum(y == 1L) / sum(y == 1L)), check.names = FALSE)
}
build_scores <- function(matrix_path, material_name, locked_all) {
  mat <- readRDS(matrix_path)
  stopifnot(is.matrix(mat), nrow(mat) == 20000L)
  sample_id <- colnames(mat)
  positive_label <- as.integer(grepl("patient", sample_id))
  conventional <- colMeans(log2(mat + 1))
  locked_sub <- locked_all[locked_all$material == material_name, , drop = FALSE]
  locked <- locked_sub[match(sample_id, locked_sub$sample_id), , drop = FALSE]
  if (anyNA(locked$frozen_score) || !identical(locked$sample_id, sample_id)) stop("Locked-score sample order does not match ", material_name)
  rbind(
    data.frame(material = material_name, sample_id, positive_label, strategy = "Conventional", score = conventional, check.names = FALSE),
    data.frame(material = material_name, sample_id, positive_label, strategy = "Locked (A/B model)", score = locked$frozen_score, check.names = FALSE)
  )
}

locked_all <- read.delim(locked_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
scores <- rbind(build_scores(cf_path, "cfDNA", locked_all), build_scores(pbl_path, "PBL", locked_all))
group_id <- interaction(scores$material, scores$strategy, drop = TRUE)
scores$auc <- NA_real_
for (g in levels(group_id)) {
  idx <- which(group_id == g)
  scores$auc[idx] <- auc_rank(scores$score[idx], scores$positive_label[idx])
}
write.table(scores, out_scores, sep = "\t", row.names = FALSE, quote = FALSE)

benchmark <- read.delim(bench_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
benchmark$strategy_display <- ifelse(benchmark$strategy == "conventional", "Conventional",
  ifelse(benchmark$strategy == "locked_A_model", "Locked (A/B model)", "Adaptive (D labels)"))
check <- unique(scores[, c("material", "strategy", "auc")])
names(check)[2:3] <- c("strategy_display", "reconstructed_auc")
check <- merge(check, benchmark[benchmark$strategy_display %in% check$strategy_display, c("material", "strategy_display", "auc")], by = c("material", "strategy_display"))
if (any(abs(check$reconstructed_auc - check$auc) > 1e-10)) stop("Reconstructed AUC does not exactly reproduce the verified benchmark.")

roc <- do.call(rbind, lapply(split(scores, group_id), roc_table))
strategy_levels <- c("Conventional", "Locked (A/B model)")
roc$strategy <- factor(roc$strategy, levels = strategy_levels)
auc_labels <- unique(scores[, c("material", "strategy", "auc")])

# Display-only monotone smoothing. Empirical AUC values above remain unchanged.
smooth_one <- function(d) {
  ux <- sort(unique(d$fpr))
  uy <- vapply(ux, function(z) max(d$tpr[d$fpr == z]), numeric(1))
  grid <- seq(0, 1, length.out = 201)
  if (length(ux) >= 3L) {
    y <- splinefun(ux, uy, method = "hyman")(grid)
  } else {
    y <- approx(ux, uy, xout = grid, rule = 2)$y
  }
  y <- cummax(pmin(pmax(y, 0), 1))
  data.frame(material = d$material[1], strategy = d$strategy[1], fpr = grid, tpr = y)
}
roc_smooth <- do.call(rbind, lapply(split(roc, interaction(roc$material, roc$strategy, drop = TRUE)), smooth_one))
roc_smooth$strategy <- factor(roc_smooth$strategy, levels = strategy_levels)
colours <- c("Conventional" = "#6C757D", "Locked (A/B model)" = "#2C638F")
theme_pub <- theme_minimal(base_family = "Arial", base_size = 16) +
  theme(text = element_text(colour = "#1C2530", family = "Arial"), axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 15, colour = "#1C2530"), axis.line = element_line(colour = "#1C2530", linewidth = 0.6),
        axis.ticks = element_line(colour = "#1C2530", linewidth = 0.4), panel.grid = element_blank(),
        panel.border = element_blank(), plot.margin = margin(10, 10, 10, 10))
make_roc_panel <- function(material_name, y_title, panel_title) {
  d <- roc_smooth[roc_smooth$material == material_name, , drop = FALSE]
  a <- auc_labels[auc_labels$material == material_name, , drop = FALSE]
  label <- sprintf("Conventional = %.3f\nLocked = %.3f", a$auc[a$strategy == "Conventional"], a$auc[a$strategy == "Locked (A/B model)"])
  ggplot(d, aes(fpr, tpr, colour = strategy)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "#9BA3A8", linewidth = 0.55) +
    geom_line(linewidth = 1.15, lineend = "round") +
    scale_colour_manual(values = colours) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    coord_fixed() + labs(title = panel_title, subtitle = label, x = "1 − Specificity", y = y_title, colour = NULL) +
    theme_pub + theme(legend.position = "bottom", plot.title = element_text(size = 19, face = "bold", hjust = 0.5, margin = margin(b = 6)),
                      plot.subtitle = element_text(size = 8.4, hjust = 0.5, lineheight = 1.0, margin = margin(b = 10)))
}
p <- make_roc_panel("cfDNA", "Sensitivity", "cfDNA (n = 50)") | make_roc_panel("PBL", "Sensitivity", "PBL (n = 49)")
p <- p + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(size = 15), legend.key.width = unit(1.7, "cm"))
ggsave(file.path(fig_dir, "Fig_6_B_v11.png"), p, width = 9.5, height = 5.7, units = "in", dpi = 600, bg = "white", device = "png")
writeLines(c("PANEL B: COHORT D STRATEGY COMPARISON", "Panel B displays empirical ROC curves, not mean AUC values.",
             "Every curve uses all scored samples in cfDNA (n=50) or PBL (n=49).", "Conventional and locked curves are deployment comparators.",
             "Adaptive uses Cohort D labels; its single benchmark AUC is retained in the table but is intentionally not promoted as a primary ROC curve.",
             "Curves use monotone smoothing for display only; all reported AUC values are unchanged empirical values from the unsmoothed sample scores.",
             "AUC values were exactly validated against D_strategy_benchmark_v1.tsv before export."),
           file.path(fig_dir, "Fig_6_B_v11_notes.txt"))
cat("FIG6_B_ROC_V2_PASS\n")
