#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."), mustWork = TRUE)
out_main <- file.path(root, "figures", "Fig_6")
out_sup  <- file.path(root, "figures", "Fig_S6")
dir.create(out_main, recursive = TRUE, showWarnings = FALSE)
dir.create(out_sup, recursive = TRUE, showWarnings = FALSE)

data_root <- file.path(root, "intermediate_data")
bench_path <- file.path(data_root, "deployment", "D_strategy_benchmark_v1.tsv")
dscore_path <- file.path(data_root, "deployment", "D_frozen_scores_v1.tsv")
cscore_path <- file.path(data_root, "deployment", "C_frozen_scores_v1.tsv")
cdf_path <- file.path(data_root, "deployment", "frozen_deployment_summary_v1.tsv")
corr_path <- file.path(data_root, "cohort_D", "derived_v1", "D_cfDNA_PBL_pairwise_correlation_v1.tsv")
fboot_path <- file.path(data_root, "cohort_F", "qc", "Cohort_F_frozen_deployment_bootstrap_v1.tsv")
fscores_path <- file.path(data_root, "cohort_F", "qc", "Cohort_F_frozen_model_scores_all_v1.tsv")
fsens_path <- file.path(data_root, "cohort_F", "qc", "Cohort_F_sensitivity_scores_v1.tsv")

required <- c(bench_path, dscore_path, cscore_path, cdf_path, corr_path, fboot_path, fscores_path, fsens_path)
if (!all(file.exists(required))) {
  stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))
}

bench <- fread(bench_path)
dscore <- fread(dscore_path)
cscore <- fread(cscore_path)
cdf <- fread(cdf_path)
corr <- fread(corr_path)
fboot <- fread(fboot_path)
fscores <- fread(fscores_path)
fsens <- fread(fsens_path)

stopifnot(all(c("material", "strategy", "auc") %in% names(bench)))
stopifnot(all(c("material", "group", "frozen_score") %in% names(dscore)))
stopifnot(all(c("material", "group", "frozen_score", "auc") %in% names(cscore)))
stopifnot(all(c("cohort", "material", "n", "auc") %in% names(cdf)))
stopifnot("spearman" %in% names(corr), nrow(corr) == 48L)
stopifnot(all(c("group", "n_total", "auc", "ci_low", "ci_high") %in% names(fboot)))
stopifnot(all(c("sample", "group", "model_score") %in% names(fscores)))
stopifnot(all(c("sample", "score_len20_500") %in% names(fsens)))

blue <- "#2C638F"
gold <- "#C98B25"
red <- "#B51E35"
grey <- "#5C6770"
light_grey <- "#E7E9EB"
dark <- "#1C2530"

theme_pub <- theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    text = element_text(colour = dark, family = "Arial"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12, colour = dark),
    axis.line = element_line(colour = dark, linewidth = 0.5),
    axis.ticks = element_line(colour = dark, linewidth = 0.4),
    panel.grid.major = element_line(colour = light_grey, linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = dark, fill = NA, linewidth = 0.6),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.margin = margin(8, 10, 8, 10)
  )

save_panel <- function(plot_obj, path, width = 5.5, height = 5.5) {
  ggsave(path, plot_obj, width = width, height = height, units = "in", dpi = 600, bg = "white", device = "png")
}

# Panel A: primary frozen cfDNA deployment in C versus D.
primary <- cdf[cohort %in% c("C", "D") & material == "cfDNA", .(cohort, n, auc)]
primary[, cohort := factor(cohort, levels = c("C", "D"))]
pA <- ggplot(primary, aes(cohort, auc)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.55, colour = grey) +
  geom_point(size = 5.5, shape = 21, fill = blue, colour = blue) +
  geom_text(aes(label = sprintf("AUC = %.3f\nn = %d", auc, n)), vjust = -1.0, size = 4.4, family = "Arial", colour = blue) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = "Deployment cohort", y = "Frozen-model AUC") +
  theme_pub

# Panel A v3: every frozen score in C and D, with cohort-level AUC retained.
score_cd <- rbindlist(list(
  cscore[material == "cfDNA", .(cohort, group, frozen_score, auc)],
  dscore[material == "cfDNA", .(cohort, group, frozen_score, auc)]
))
score_cd[, cohort := factor(cohort, levels = c("C", "D"))]
score_cd[, group_display := fifelse(group %in% c("healthy", "NCC"), "Healthy", group)]
score_cd[, group_display := factor(group_display, levels = c("Healthy", "SCLC", "HNSCC"))]
ann_cd <- score_cd[, .(n = .N, auc = unique(auc)[1]), by = cohort]
ann_cd[, label := sprintf("AUC = %.3f\nn = %d", auc, n)]
ann_cd[, x := 1.5]
pA_v3 <- ggplot(score_cd, aes(group_display, frozen_score, fill = group_display, colour = group_display)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.38, linewidth = 0.55) +
  geom_jitter(width = 0.10, height = 0, size = 1.65, alpha = 0.68) +
  geom_text(data = ann_cd, aes(x = x, y = 1.04, label = label), inherit.aes = FALSE, colour = dark, family = "Arial", size = 3.8, vjust = 1) +
  facet_grid(. ~ cohort, scales = "free_x") +
  scale_fill_manual(values = c(Healthy = "#B9C0C5", SCLC = blue, HNSCC = blue), drop = FALSE) +
  scale_colour_manual(values = c(Healthy = "#6C757D", SCLC = blue, HNSCC = blue), drop = FALSE) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = NULL, y = "Frozen-model score") +
  theme_pub + theme(legend.position = "none", strip.text = element_text(size = 13, face = "bold"))

# Panel A clean version: retain every sample point and distribution, remove only
# the cohort-level AUC/n annotation for a less crowded standalone panel.
pA_v4 <- pA_v3
text_layers <- vapply(pA_v4$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1))
pA_v4$layers <- pA_v4$layers[!text_layers]

# Panel B: D conventional/adaptive/locked benchmark; adaptive is explicitly a comparator.
bench[, strategy := factor(strategy, levels = c("conventional", "adaptive_5fold", "locked_A_model"))]
bench[, material := factor(material, levels = c("cfDNA", "PBL"))]
bench_lab <- c(conventional = "Conventional", adaptive_5fold = "Adaptive\n(D labels)", locked_A_model = "Locked\n(A/B model)")
pB <- ggplot(bench, aes(strategy, auc, colour = material)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.55, colour = grey) +
  geom_point(position = position_dodge(width = 0.35), size = 4.8) +
  geom_text(aes(label = sprintf("%.3f", auc)), position = position_dodge(width = 0.35), vjust = -1.0, size = 3.7, family = "Arial", show.legend = FALSE) +
  scale_x_discrete(labels = bench_lab) +
  scale_colour_manual(values = c(cfDNA = blue, PBL = gold)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.04))) +
  labs(x = "Deployment strategy", y = "AUC", colour = NULL) +
  theme_pub

# Panel C: D sample-level frozen score distributions.
dscore[, material := factor(material, levels = c("cfDNA", "PBL"))]
dscore[, group := factor(group, levels = c("healthy", "HNSCC"))]
pC <- ggplot(dscore, aes(group, frozen_score)) +
  geom_boxplot(aes(fill = material), outlier.shape = NA, width = 0.55, alpha = 0.32, colour = dark, linewidth = 0.55) +
  geom_jitter(aes(colour = material), width = 0.10, height = 0, size = 1.7, alpha = 0.65) +
  facet_grid(. ~ material, scales = "free_x", labeller = as_labeller(c(cfDNA = "cfDNA", PBL = "PBL"))) +
  scale_fill_manual(values = c(cfDNA = blue, PBL = gold)) +
  scale_colour_manual(values = c(cfDNA = blue, PBL = gold)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = NULL, y = "Frozen-model score") +
  theme_pub + theme(legend.position = "none", strip.text = element_text(size = 13, face = "bold"))

# Panel D: all paired rows present in the current source table; no binning.
corr[, x := 1]
med_rho <- median(corr$spearman, na.rm = TRUE)
pD <- ggplot(corr, aes(x, spearman)) +
  geom_jitter(width = 0.12, height = 0, size = 3.0, shape = 21, fill = "white", colour = blue, stroke = 1.1) +
  geom_hline(yintercept = med_rho, colour = red, linewidth = 1.0) +
  annotate("text", x = 1.08, y = 0.856, label = sprintf("median rho = %.3f", med_rho), hjust = 1, vjust = 1, colour = red, size = 4.4, family = "Arial") +
  scale_x_continuous(limits = c(0.72, 1.28), breaks = NULL) +
  scale_y_continuous(limits = c(0.65, 0.86), breaks = seq(0.65, 0.85, 0.05), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = NULL, y = "Paired cfDNA–PBL Spearman correlation") +
  theme_pub + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Panel D v3: ranked lollipop display preserves every paired observation and exposes low-correlation pairs.
corr_rank <- copy(corr)[order(spearman)]
corr_rank[, pair_rank := seq_len(.N)]
pD_v3 <- ggplot(corr_rank, aes(spearman, pair_rank)) +
  geom_vline(xintercept = med_rho, colour = red, linewidth = 1.0) +
  geom_segment(aes(x = 0.65, xend = spearman, y = pair_rank, yend = pair_rank), colour = "#B8C1C8", linewidth = 0.55) +
  geom_point(shape = 21, size = 3.1, stroke = 1.0, fill = "white", colour = blue) +
  annotate("text", x = med_rho + 0.002, y = 47.5, label = sprintf("median rho = %.3f", med_rho), hjust = 0, vjust = 1, colour = red, size = 4.1, family = "Arial") +
  scale_x_continuous(limits = c(0.65, 0.86), breaks = seq(0.65, 0.85, 0.05), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(0.5, 48.5), breaks = c(1, 12, 24, 36, 48), expand = c(0, 0)) +
  labs(x = "Paired cfDNA–PBL Spearman correlation", y = "Pair rank (low to high)") +
  theme_pub + theme(panel.grid.major.y = element_line(colour = light_grey, linewidth = 0.25))

# Panel E: F cancer-stratified frozen deployment with bootstrap intervals.
fboot[, group_label := paste0(group, " (n=", n_total, ")")]
fboot[, group_label := factor(group_label, levels = group_label[order(auc)])]
pE <- ggplot(fboot, aes(auc, group_label)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.55, colour = grey) +
  geom_segment(aes(x = ci_low, xend = ci_high, y = group_label, yend = group_label), linewidth = 1.1, colour = blue) +
  geom_point(shape = 21, size = 4.2, stroke = 1.1, fill = "white", colour = blue) +
  scale_x_continuous(limits = c(0.2, 1.02), breaks = seq(0.2, 1, 0.2), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Frozen-model AUC versus Healthy", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 10), panel.grid.major.y = element_line(colour = light_grey, linewidth = 0.3))

# Panel S6-A: F integrity and eligibility audit. This is a scope/QC panel, not a performance estimate.
flow_df <- data.frame(
  x = c(1.2, 4.0, 6.8, 9.5),
  label = c("301\nMeDIP BED files\nprocessed", "283\nbaseline proxy\nunits", "18\ntrailing b\nreplicates excluded", "20,000\nfeatures in exact\nfrozen order")
)
pS_A <- ggplot() +
  annotate("rect", xmin = 0.45, xmax = 2.0, ymin = 4.0, ymax = 6.0, fill = "#EAF1F6", colour = blue, linewidth = 0.8) +
  annotate("rect", xmin = 3.25, xmax = 4.75, ymin = 4.0, ymax = 6.0, fill = "#EAF1F6", colour = blue, linewidth = 0.8) +
  annotate("rect", xmin = 6.05, xmax = 7.55, ymin = 4.0, ymax = 6.0, fill = "#FBF0D7", colour = gold, linewidth = 0.8) +
  annotate("rect", xmin = 8.75, xmax = 10.25, ymin = 4.0, ymax = 6.0, fill = "#EAF1F6", colour = blue, linewidth = 0.8) +
  annotate("segment", x = 2.0, xend = 3.25, y = 5.0, yend = 5.0, linewidth = 0.8, arrow = arrow(length = unit(0.18, "cm"), type = "closed"), colour = grey) +
  annotate("segment", x = 4.75, xend = 6.05, y = 5.0, yend = 5.0, linewidth = 0.8, arrow = arrow(length = unit(0.18, "cm"), type = "closed"), colour = grey) +
  annotate("segment", x = 7.55, xend = 8.75, y = 5.0, yend = 5.0, linewidth = 0.8, arrow = arrow(length = unit(0.18, "cm"), type = "closed"), colour = grey) +
  geom_text(data = flow_df, aes(x, 5, label = label), family = "Arial", size = 4.0, lineheight = 0.95, colour = dark) +
  annotate("text", x = 5.35, y = 2.8, label = "No F-specific feature selection, scaling recalibration, model fitting, or threshold optimization", family = "Arial", size = 4.0, colour = dark) +
  coord_cartesian(xlim = c(0, 10.7), ylim = c(2.1, 7.1), expand = FALSE) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

# Panel S6-B: fragment-length sensitivity pilot, limited to the six audited pilot records.
fsens[, sample := factor(sample, levels = sample)]
pS_B <- ggplot(fsens, aes(sample, score_len20_500)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.55, colour = grey) +
  geom_point(size = 3.4, shape = 21, fill = blue, colour = blue) +
  geom_text(aes(label = sprintf("%.3f", score_len20_500)), vjust = -1.0, size = 3.2, family = "Arial", colour = blue) +
  scale_y_continuous(limits = c(0.45, 1.04), breaks = seq(0.5, 1, 0.1), expand = expansion(mult = c(0.01, 0.03))) +
  labs(x = "Pilot GEO record", y = "Frozen score under 20–500 bp sensitivity policy") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# Panel S6-C: binned view of the same D correlation table for a supplementary shape check.
pS_C <- ggplot(corr, aes(spearman)) +
  geom_histogram(binwidth = 0.01, boundary = 0.65, fill = blue, colour = "white", linewidth = 0.35) +
  geom_vline(xintercept = med_rho, colour = red, linewidth = 0.9) +
  annotate("text", x = med_rho + 0.002, y = Inf, label = sprintf("median rho = %.3f", med_rho), hjust = 0, vjust = 1.2, colour = red, size = 3.6, family = "Arial") +
  scale_x_continuous(limits = c(0.65, 0.86), breaks = seq(0.65, 0.85, 0.05), expand = expansion(mult = c(0.01, 0.01))) +
  labs(x = "Paired cfDNA–PBL Spearman correlation", y = "Pairs") +
  theme_pub

# Panel S6-D: F sample-level frozen-score distributions, grouped by cancer type.
fmed <- fscores[, .(med = median(model_score, na.rm = TRUE)), by = group][order(med)]
fscores[, group := factor(group, levels = fmed$group)]
pS_D <- ggplot(fscores, aes(model_score, group)) +
  geom_boxplot(outlier.shape = NA, fill = "#D9E5EE", colour = blue, linewidth = 0.55, width = 0.6) +
  geom_jitter(height = 0.10, width = 0, size = 0.9, alpha = 0.45, colour = blue) +
  scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Frozen C model score", y = NULL) +
  theme_pub + theme(axis.text.y = element_text(size = 9))

save_panel(pA, file.path(out_main, "Fig_6_A_v2.png"), 5.5, 5.5)
save_panel(pA_v3, file.path(out_main, "Fig_6_A_v3.png"), 5.5, 5.5)
save_panel(pA_v4, file.path(out_main, "Fig_6_A_v4.png"), 5.5, 5.5)
save_panel(pB, file.path(out_main, "Fig_6_B_v2.png"), 5.5, 5.5)
save_panel(pC, file.path(out_main, "Fig_6_C_v2.png"), 5.5, 5.5)
save_panel(pD, file.path(out_main, "Fig_6_D_v2.png"), 5.5, 5.5)
save_panel(pD_v3, file.path(out_main, "Fig_6_D_v3.png"), 5.5, 5.5)
save_panel(pE, file.path(out_main, "Fig_6_E_v2.png"), 8.8, 5.5)

main_fig_v3 <- (pA_v3 | pB) / (pC | pD) / pE +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 6. Locked deployment, strategy dependence and external portability across cohorts C–F",
    subtitle = "C and D provide independent cfDNA deployment and paired-background stress tests; F provides cross-center, cross-cancer deployment in the same frozen representation.",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 18, colour = dark),
      plot.subtitle = element_text(family = "Arial", size = 12, colour = grey),
      plot.tag = element_text(family = "Arial", face = "bold", size = 17, colour = dark),
      plot.margin = margin(8, 8, 8, 8)
    )
  )
ggsave(file.path(out_main, "Fig_6_natural_v3.png"), main_fig_v3, width = 15, height = 17, units = "in", dpi = 600, bg = "white", device = "png")

main_fig_v4 <- (pA_v3 | pB) / (pC | pD_v3) / pE +
  plot_annotation(
    tag_levels = "A",
    title = "Figure 6. Locked deployment, strategy dependence and external portability across cohorts C–F",
    subtitle = "C and D provide independent cfDNA deployment and paired-background stress tests; F provides cross-center, cross-cancer deployment in the same frozen representation.",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 18, colour = dark),
      plot.subtitle = element_text(family = "Arial", size = 12, colour = grey),
      plot.tag = element_text(family = "Arial", face = "bold", size = 17, colour = dark),
      plot.margin = margin(8, 8, 8, 8)
    )
  )
ggsave(file.path(out_main, "Fig_6_natural_v4.png"), main_fig_v4, width = 15, height = 17, units = "in", dpi = 600, bg = "white", device = "png")

save_panel(pS_A, file.path(out_sup, "Fig_S6_A_v1.png"), 7.2, 4.5)
save_panel(pS_B, file.path(out_sup, "Fig_S6_B_v1.png"), 6.5, 5.0)
save_panel(pS_C, file.path(out_sup, "Fig_S6_C_v1.png"), 6.5, 5.0)
save_panel(pS_D, file.path(out_sup, "Fig_S6_D_v1.png"), 8.0, 6.0)

sup_fig <- (pS_A | pS_B) / (pS_C | pS_D) +
  plot_annotation(
    tag_levels = "A",
    title = "Supplementary Figure S6. Cohort D/F integrity, sensitivity and score-distribution details",
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 17, colour = dark),
      plot.tag = element_text(family = "Arial", face = "bold", size = 16, colour = dark),
      plot.margin = margin(8, 8, 8, 8)
    )
  )
ggsave(file.path(out_sup, "Fig_S6_natural_v1.png"), sup_fig, width = 15, height = 13, units = "in", dpi = 600, bg = "white", device = "png")

main_legend <- c(
  "FIGURE TITLE",
  "Locked deployment, strategy dependence and external portability across cohorts C–F.",
  "",
  "WHAT THE FIGURE SHOWS",
  "A, every frozen cfDNA score in the independent C and D deployment cohorts, stratified by Healthy versus SCLC or HNSCC; AUC and total n are shown per cohort. B, Cohort D conventional, D-label-adaptive comparator and A/B-locked strategy benchmark for cfDNA and paired PBL. C, sample-level frozen-score distributions in D, stratified by material and healthy/HNSCC group. D, all 48 rows in the current D paired cfDNA–PBL correlation source table, with the median Spearman correlation marked in red. E, Cohort F cancer-stratified frozen-model AUC estimates versus Healthy with stratified bootstrap 95% intervals.",
  "",
  "KEY RESULTS",
  "C cfDNA AUC = 0.864 (n = 94); D cfDNA AUC = 0.663 (n = 50). In D, locked cfDNA AUC = 0.663 versus conventional = 0.493 and D-adaptive = 0.400; locked PBL AUC = 0.507 and is treated as a secondary failure boundary. The D paired-correlation table contains 48 rows and has median rho = 0.796. F contains 301 processed BED records, 283 title-derived baseline proxy units after excluding 18 trailing-b replicates, and 15 cancer-group versus Healthy estimates.",
  "",
  "INTERPRETATION",
  "The figure separates three questions: primary frozen cfDNA deployment, strategy dependence under the HNSCC stress test and paired leukocyte-background concordance. Cohort F is shown as cancer-stratified portability evidence rather than a pooled performance estimate. Near-chance F groups remain visible as deployment boundaries and are not removed.",
  "",
  "DATA AND METHODS",
  paste0("Inputs: ", bench_path, "; ", dscore_path, "; ", cdf_path, "; ", corr_path, "; ", fboot_path, "; ", fscores_path, "; ", fsens_path, ". A/B locked features, order, scaling, model coefficients and threshold were held fixed for C, D and F. No F-specific feature selection, scaling recalibration, model fitting or threshold optimization was performed. PBL_control_20 was excluded from D PBL analysis after the repository duplication audit."),
  "",
  "EVIDENCE LEVEL AND LIMITATIONS",
  "C and D are independent deployment cohorts with different sample structures. F subject units are title-derived baseline proxies because GEO metadata does not provide a complete explicit subject/timepoint key; this limitation must remain explicit. The D correlation source table contains 48 rows despite an earlier 49-pair summary. Adaptive results use D labels and are comparator analyses, not locked deployment claims. Intervals are stratified bootstrap 95% intervals for F cancer-group estimates.",
  "",
  "COLOR KEY",
  "Blue = cfDNA or frozen-model deployment; gold = PBL in the D strategy benchmark; red = D paired-correlation median; dashed grey = chance AUC = 0.5."
)
writeLines(main_legend, file.path(out_main, "Fig_6_title_legend_v3.txt"))
main_legend_v4 <- sub("all 48 rows in the current D paired cfDNA–PBL correlation source table, with the median Spearman correlation marked in red", "all 48 rows in the current D paired cfDNA–PBL correlation source table, ranked from low to high; a red vertical line marks the median Spearman correlation", main_legend, fixed = TRUE)
writeLines(main_legend_v4, file.path(out_main, "Fig_6_title_legend_v4.txt"))

sup_legend <- c(
  "FIGURE TITLE",
  "Cohort D/F integrity, sensitivity and score-distribution details.",
  "",
  "WHAT THE FIGURE SHOWS",
  "A, Cohort F processing and eligibility audit. B, six-record fragment-length sensitivity pilot under the 20–500 bp policy. C, binned display of all rows in the D paired-correlation table. D, sample-level frozen-score distributions across Healthy and F cancer groups.",
  "",
  "KEY RESULTS",
  "F processing retained 301 BED records for representation construction, collapsed them to 283 title-derived baseline proxy units, excluded 18 trailing-b replicates and generated exact 20,000-feature representations. The six-record sensitivity pilot is descriptive only. The D histogram is a display alternative to the point distribution in the main figure.",
  "",
  "INTERPRETATION",
  "These panels document input integrity, sensitivity context and the score distributions underlying the primary deployment summaries; they do not add cohort-adaptive analyses.",
  "",
  "DATA AND METHODS",
  paste0("Inputs: ", fsens_path, "; ", corr_path, "; ", fscores_path, ". The frozen 20–1000 bp policy and exact C feature order were used for F representation construction. The fragment-length panel is restricted to the six records in the sensitivity table."),
  "",
  "EVIDENCE LEVEL AND LIMITATIONS",
  "F subject identifiers are title-derived proxy units, not confirmed independent subjects. The sensitivity panel is a six-record pilot and must not be interpreted as a cohort-wide robustness estimate. D uses the current 48-row correlation source table.",
  "",
  "COLOR KEY",
  "Blue = frozen representation/model outputs; red = paired-correlation median; dashed grey = chance AUC = 0.5."
)
writeLines(sup_legend, file.path(out_sup, "Fig_S6_title_legend_v1.txt"))

cat("FIGURE6_EXPANDED_PASS\n")
cat("main:", file.path(out_main, "Fig_6_natural_v4.png"), "\n")
cat("supplement:", file.path(out_sup, "Fig_S6_natural_v1.png"), "\n")
