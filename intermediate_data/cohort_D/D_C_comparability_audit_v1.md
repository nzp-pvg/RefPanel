# Cohort C/D comparability audit (v1)

## Checks completed

| Check | Finding | Consequence |
|---|---|---|
| D numeric field | BED column 4 is integer-like raw signal/count values; observed locked-matrix range is 0–1792 for cfDNA and 0–1709 for PBL. | D summaries must use the same `log2(value+1)` transform documented for C; raw D values must not be compared directly with C log2 matrices. |
| C transform | The public C pipeline explicitly reads raw count matrices and applies `log2(count+1)`. | The D transfer summaries used the same transform, so the current descriptive signal scale is aligned at the transform level. |
| Library-size normalization | No TMM/CPM/library-size normalization is applied in the public C preprocessing block. | This is a remaining cross-cohort technical limitation; sample-level total signal and nonzero fraction must remain QC covariates, not biological conclusions. |
| Coordinates | A/D intervals are 1-based inclusive 300-bp intervals (`1–300`, `301–600`); locked features are 300 bp. | Coordinate matching is internally consistent, but files must not be passed to tools expecting standard 0-based half-open BED without conversion. |
| Feature provenance | D matrices use the A/B-only 705,939-feature lock; no D variance ranking or threshold tuning was used. | No target-cohort feature-selection leakage in the current D extraction. |
| PBL exclusion | `PBL_control_20` is excluded from the 49-sample secondary matrix and retained for audit. | Secondary PBL result is not contaminated by the known duplicate file. |
| C/D disease context | C is SCLC/NCC; D is HNSCC/healthy. | D is a cross-disease stress test, not a same-disease replication cohort. |

## Go/no-go decision

The current D result is interpretable as a frozen cross-disease transfer stress test at the transform level, but not as a definitive test of clinical generalization because the cohorts differ in disease, acquisition context and likely signal distribution. The primary cfDNA result shows no useful separation in this frozen representation (Cohen d −0.082; PC1 descriptive AUC 0.463). Therefore the original plan to claim performance improvement from adding D is a no-go. The defensible endpoint is a negative/neutral external-transfer result plus a reproducibility demonstration.

Further work should be limited to pre-specified QC and paired-shift summaries. Do not reselect features, refit thresholds, or harmonize distributions using D labels. A positive claim would require a new development cohort or an independently pre-specified model trained before D was opened.
