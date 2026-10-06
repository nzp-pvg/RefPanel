# Cohort F deployment gate v2

| criterion | evidence | status |
|---|---|---|
| all public MeDIP BED files | 301/301 processed | PASS |
| malformed BED rows | 0 | PASS |
| exact frozen C feature space | 301/301 matrices, 20,000 features | PASS |
| baseline/replicate adjudication | 283 baseline proxy samples; 18 `b` replicates excluded | PASS WITH PROXY LIMITATION |
| frozen C scores | 301/301 finite scores | PASS |
| bootstrap uncertainty | 2,000 stratified resamples per cancer type | PASS |
| external feature/model adaptation | none | PASS |

The cohort is reportable as a frozen cross-center, cross-cancer deployment cohort. Subject identifiers are adjudicated title-derived proxy units; this limitation must remain explicit. No feature selection, model fitting, scaling recalibration, or threshold optimization was performed in F.
