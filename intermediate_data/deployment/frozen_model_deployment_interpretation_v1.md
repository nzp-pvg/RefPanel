# Frozen model deployment checkpoint (v1)

The model was trained only on Cohort A HCC versus CTL using 20,000 A-derived features, L1 logistic regression (`cv.glmnet`, alpha=1, 5-fold CV, lambda.1se=0.0931416, internal standardization). Cohort C and Cohort D were scored without feature retuning.

| Deployment | Samples | Frozen-model AUC | Interpretation |
|---|---:|---:|---|
| C cfDNA, SCLC vs NCC | 74 + 20 | 0.8642 | Cross-disease/domain deployment signal remains, but this is not the archived within-C CV estimate. |
| D cfDNA, HNSCC vs healthy | 30 + 20 | 0.6633 | Weak-to-moderate exploratory portability signal; below a confirmatory threshold and not evidence of improvement. |
| D PBL, HNSCC vs healthy | 30 + 19 | 0.5070 | No frozen-model discrimination. |

The C/D scores are a deployment test of a model trained in A, not a newly optimized classifier. D therefore supports a domain-stress-test interpretation: some score portability is present in cfDNA, but it is attenuated relative to C and absent in PBL. The next benchmark should compare this locked result with conventional and adaptive comparators, while keeping the locked result primary and declaring a deployment-failure flag for D PBL.
