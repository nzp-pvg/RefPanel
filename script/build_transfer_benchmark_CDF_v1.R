#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "project_paths.R"))
outdir <- file.path(data_root, 'deployment')
d <- fread(file.path(outdir,'D_strategy_benchmark_v1.tsv'))
f <- fread(file.path(data_root, 'cohort_F', 'qc', 'Cohort_F_frozen_model_scores_all_v1.tsv'))
fm <- fread(file.path(data_root, 'cohort_F', 'derived', 'GSE243474_locked_manifest_v1.tsv'))
f <- merge(f, fm[,.(sample=geo_accession,primary_eligible)], by='sample')
f <- f[primary_eligible=='YES',]
summ <- f[, .(n=.N, mean_score=mean(model_score), median_score=median(model_score), mean_nonzero=mean(nonzero_fraction)), by=group]
fwrite(d,file.path(outdir,'C_D_F_strategy_benchmark_v1.tsv'),sep='\t')
fwrite(summ,file.path(data_root, 'cohort_F', 'qc', 'Cohort_F_group_deployment_summary_v1.tsv'),sep='\t')
writeLines(c('# C/D/F benchmark and deployment summary v1','D strategy benchmark is retained as the explicit adaptive-vs-locked comparator. Cohort F group summaries are descriptive locked deployment outputs; no F-adaptive model was fit.','The F manifest and bootstrap table define eligibility and uncertainty.'), file.path(outdir, 'C_D_F_benchmark_notes_v1.md'))
