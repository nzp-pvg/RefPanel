#!/usr/bin/env Rscript

# Cohort D cross-material heatmap in the exact frozen 20,000-feature model space.
# This script intentionally performs no model fitting, feature selection, or scaling update.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "..", ".."), mustWork = TRUE)
data_root <- file.path(root, "intermediate_data")
derived_dir <- file.path(data_root, "cohort_D", "derived_v1")
fig_dir <- file.path(root, "figures", "Fig_6")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cf_rds <- file.path(derived_dir, "D_cfDNA_model20k_v1.rds")
pbl_rds <- file.path(derived_dir, "D_PBL_model20k_v1.rds")
pairs_path <- file.path(derived_dir, "D_cfDNA_PBL_pairwise_correlation_v1.tsv")
out_matrix <- file.path(derived_dir, "D_cfDNA_PBL_cross_correlation_model20k_v1.tsv")

required <- c(cf_rds, pbl_rds, pairs_path)
if (!all(file.exists(required))) stop("Missing required input(s): ", paste(required[!file.exists(required)], collapse = ", "))

# Frozen input contract: exact model matrices and the already-audited 48 matched pairs only.
cf <- readRDS(cf_rds)
pbl <- readRDS(pbl_rds)
pairs <- fread(pairs_path)
stopifnot(is.matrix(cf), is.matrix(pbl), nrow(cf) == 20000L, nrow(pbl) == 20000L)
stopifnot(identical(rownames(cf), rownames(pbl)))
stopifnot(identical(names(pairs), c("subject", "cfDNA", "PBL", "spearman")), nrow(pairs) == 48L)
stopifnot(!anyDuplicated(pairs$subject), !anyDuplicated(pairs$cfDNA), !anyDuplicated(pairs$PBL))
stopifnot(all(pairs$cfDNA %in% colnames(cf)), all(pairs$PBL %in% colnames(pbl)))

# Pair order is fixed from the manifest-derived pairwise source, so the diagonal is the matched individual.
cf_pair <- cf[, pairs$cfDNA, drop = FALSE]
pbl_pair <- pbl[, pairs$PBL, drop = FALSE]
rho <- cor(cf_pair, pbl_pair, method = "spearman", use = "pairwise.complete.obs")
stopifnot(identical(dim(rho), c(48L, 48L)), all(is.finite(rho)))

diag_model20k <- diag(rho)
diag_check <- data.table(subject = pairs$subject, rho_model20k = diag_model20k, rho_full_locked = pairs$spearman)
stopifnot(nrow(diag_check) == 48L, !anyNA(diag_check$rho_model20k))

rho_long <- as.data.table(as.table(rho))
setnames(rho_long, c("cfDNA", "PBL", "rho_model20k"))
rho_long[, cf_pair_index := match(cfDNA, pairs$cfDNA)]
rho_long[, pbl_pair_index := match(PBL, pairs$PBL)]
rho_long[, matched_pair := cf_pair_index == pbl_pair_index]
setorder(rho_long, cf_pair_index, pbl_pair_index)
fwrite(rho_long, out_matrix, sep = "\t")

blue <- "#2C638F"
red <- "#B51E35"
dark <- "#1C2530"
light_grey <- "#E7E9EB"
theme_pub <- theme_minimal(base_family = "Arial", base_size = 14) +
  theme(
    text = element_text(colour = dark, family = "Arial"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 10, colour = dark),
    axis.line = element_line(colour = dark, linewidth = 0.5),
    axis.ticks = element_line(colour = dark, linewidth = 0.4),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = dark, fill = NA, linewidth = 0.6),
    plot.margin = margin(8, 10, 8, 10)
  )

# A full 48 x 48 matrix has no readable per-cell text at journal panel scale;
# the diagonal is outlined rather than relying only on its positional interpretation.
heat <- ggplot(rho_long, aes(pbl_pair_index, cf_pair_index, fill = rho_model20k)) +
  geom_tile(colour = "white", linewidth = 0.08) +
  geom_tile(data = rho_long[matched_pair == TRUE], fill = NA, colour = red, linewidth = 0.42) +
  scale_fill_gradientn(
    colours = c("#0B3C5D", "#2F6F94", "#9BA99B", "#FFD36F", "#E64B6A"),
    values = scales::rescale(c(0.35, 0.50, 0.65, 0.77, 0.90)),
    limits = c(0.35, 0.90), oob = scales::squish, name = expression(rho)
  ) +
  scale_x_continuous(breaks = c(1, 12, 24, 36, 48), expand = c(0, 0)) +
  scale_y_reverse(breaks = c(1, 12, 24, 36, 48), expand = c(0, 0)) +
  coord_fixed() +
  labs(x = "PBL pair index", y = "cfDNA pair index") +
  guides(fill = guide_colorbar(barheight = unit(4.2, "cm"), barwidth = unit(0.55, "cm"))) +
  theme_pub + theme(legend.position = "right")

pair_cmp <- rbindlist(list(
  data.table(type = "Matched diagonal", rho = diag_model20k),
  data.table(type = "All non-matched cells", rho = rho[row(rho) != col(rho)])
))
pair_cmp[, type := factor(type, levels = c("All non-matched cells", "Matched diagonal"))]
cmp <- ggplot(pair_cmp, aes(type, rho, fill = type)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.50, linewidth = 0.55) +
  geom_jitter(width = 0.10, height = 0, size = 0.55, alpha = 0.15, colour = blue) +
  scale_fill_manual(values = c("All non-matched cells" = "#B9C0C5", "Matched diagonal" = blue)) +
  scale_y_continuous(limits = c(0.35, 0.90), breaks = seq(0.4, 0.9, 0.1)) +
  labs(x = NULL, y = expression("Cross-material Spearman " * rho)) +
  theme_pub + theme(legend.position = "none", axis.text.x = element_text(size = 10))

save_png <- function(plot_obj, path, width, height) {
  ggsave(path, plot_obj, width = width, height = height, units = "in", dpi = 600, bg = "white", device = "png")
}

save_png(heat, file.path(fig_dir, "Fig_6_D_v5.png"), 5.5, 5.5)
save_png(cmp, file.path(fig_dir, "Fig_6_D_inset_v1.png"), 4.2, 5.5)

summary_lines <- c(
  "COHORT D CROSS-MATERIAL HEATMAP",
  "Frozen input space: exact 20,000 model features; no D-specific feature selection, scaling, or model fitting.",
  sprintf("Matrix: %d cfDNA–PBL cross-material correlations (%d matched diagonal pairs and %d non-matched cells).", nrow(rho_long), nrow(pairs), nrow(rho_long) - nrow(pairs)),
  sprintf("Median matched diagonal rho in 20,000-feature space: %.3f.", median(diag_model20k)),
  sprintf("Median non-matched rho in 20,000-feature space: %.3f.", median(rho[row(rho) != col(rho)])),
  sprintf("Correlation between matched 20,000-feature rho and full-locked rho: %.3f.", cor(diag_check$rho_model20k, diag_check$rho_full_locked, method = "spearman")),
  "Rows and columns follow the same 48-pair order; red outlined cells are matched cfDNA–PBL pairs.",
  "PBL_control_20 is absent because it was excluded by the repository duplication audit."
)
writeLines(summary_lines, file.path(fig_dir, "Fig_6_D_cross_material_heatmap_v1.txt"))

cat("CROSS_MATERIAL_HEATMAP_PASS\n")
cat("matrix:", out_matrix, "\n")
cat("panel:", file.path(fig_dir, "Fig_6_D_v5.png"), "\n")
