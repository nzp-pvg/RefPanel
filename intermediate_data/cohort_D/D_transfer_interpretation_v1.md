# Cohort D frozen-transfer checkpoint (v1)

## Status

The A/B-only 300-bp lock contains 705,939 features. Cohort D matrices were extracted without feature ranking or threshold tuning: primary cfDNA is 30 HNSCC plus 20 healthy controls; secondary PBL is 30 HNSCC plus 19 healthy controls, with `PBL_control_20` excluded because it is byte-for-byte identical to `FaDu_MeDIP.bed`.

## Results

| Material | n HNSCC | n healthy | Mean log2(signal) difference (HNSCC−healthy) | Cohen d | Wilcoxon p | PC1 descriptive AUC |
|---|---:|---:|---:|---:|---:|---:|
| cfDNA primary | 30 | 20 | −0.0291 | −0.082 | 0.6703 | 0.4633 |
| PBL secondary | 30 | 19 | 0.1417 | 0.411 | 0.5180 | 0.4456 |

Median cfDNA–PBL paired Spearman correlation across the 49 available same-number pairs is 0.7963. PCA and PC1 AUC are descriptive summaries computed in the frozen feature space; no claim of classifier performance or clinical validation is made.

## Interpretation and decision gate

The primary cfDNA transfer does not currently show a useful global separation signal: the mean effect is near zero and PC1 is close to chance. The secondary PBL analysis has a modest effect-size direction, but its uncertainty is large and the PC1 direction is not discriminatory. These data do not yet confirm an improvement potential. A limited, defensible next step is subgroup/paired-shift QC using the same locked features and pre-specified metadata only; any model fitting or threshold estimation would require a separately declared development/training design and must not be called external validation.

## Inputs and limitations

Inputs are `AB_locked_features_v1.bed`, the D sample manifest, and the deposited 1-based inclusive interval BED files. Values are the fourth BED column, transformed only as `log2(value+1)` for summaries/PCA. The analysis is limited by 50 cfDNA and 49 PBL samples, lack of a D-specific clinical outcome, possible platform/protocol differences, and the fact that the current A/B lock is reconstructed from archived candidate features and B IC matrices. The excluded PBL file remains retained for audit.
