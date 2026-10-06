# Cohort F frozen deployment analysis summary

## Data readiness

All 301 public MeDIP fragment BED files were downloaded and processed. The final deduplicated BED QC table contains 301 files. No malformed BED rows were detected. Coordinates and fragment lengths were screened with the frozen 20–1000 bp policy; flagged records were excluded from representation construction.

## Subject-level lock

The 301 records were adjudicated using the sample-title replicate convention. The 18 records with a trailing `b` replicate suffix were excluded, leaving 283 adjudicated baseline proxy samples. No Cohort F feature, scaling parameter, model parameter, or threshold was estimated from these data.

## Frozen representation and deployment

Every record has a 20,000-feature count file in the exact frozen C feature order. Frozen C scores are in `Cohort_F_frozen_model_scores_all_v1.tsv`; stratified bootstrap AUC estimates (2,000 resamples per cancer type versus Healthy) are in `Cohort_F_frozen_deployment_bootstrap_v1.tsv`.

## Interpretation boundary

These are locked deployment results. They support cross-center and cross-cancer portability assessment, not cohort-adaptive biomarker discovery. Cancer-specific estimates must be reported with their sample sizes and confidence intervals; pooled performance is not a primary endpoint.
