# Manuscript reconstruction plan for BMC Bioinformatics

## Central claim

An A/B-calibrated, input-free cfMeDIP representation can be deployed without external feature retuning, but portability is domain-dependent and should be governed by explicit deployment QC.

## Main evidence sequence

1. Cohort A defines stable genomic representation and candidate feature universe.
2. Cohort B uses paired IP–IC data for calibration only; all model and transformation settings are frozen.
3. Cohort C is the primary locked transfer evaluation.
4. Cohort D is an independent HNSCC cfDNA–PBL stress test; the PBL duplicate is excluded transparently.
5. Cohort F is the multi-cancer, cross-center frozen deployment cohort.
6. Conventional, adaptive and locked strategies are benchmarked; adaptive results are labeled exploratory.
7. QC gates define reportable versus untrusted deployment outputs.

## Figure/table package

- Figure 1: A/B development and freeze boundary.
- Figure 2: exact representation and reproducibility workflow.
- Figure 3: C/D/F locked deployment design and sample flow.
- Figure 4: cancer-stratified deployment performance and failure boundaries.
- Supplementary Table 1: frozen manifest and artifact checksums.
- Supplementary Table 2: benchmark and ablation results.
- Supplementary Table 3: D/F QC flags, replicate adjudication and exclusions.

## Prohibited claims

Do not describe D/F as adaptive discovery cohorts, do not pool cancer types as a primary endpoint, and do not interpret low-performing domains as data failures to be removed. Report them as portability boundaries.
