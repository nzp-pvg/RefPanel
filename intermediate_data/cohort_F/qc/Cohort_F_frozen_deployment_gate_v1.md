# Cohort F frozen deployment gate v1

All six pilot MeDIP fragment BED files (3 Healthy, 3 Colorectal) were converted to the same A-derived hg19 300-bp stable feature universe. The fixed representation policy retained records with `start >= 0` and `20 <= end-start <= 1000` bp; malformed and long-span records were excluded and documented.

| criterion | result | gate |
|---|---|---|
| hg19/chr coordinate compatibility | passed in all 6 files | PASS |
| exact A 300-bp feature universe | 882,903 features returned for all 6 | PASS |
| finite coverage/count vector | passed in all 6 | PASS |
| nonzero feature fraction | 0.783–0.920 | PASS for feasibility |
| independent-subject mapping | not resolved by title metadata alone | PENDING |
| fragment anomaly policy | fixed and auditable; sensitivity analysis still required | PENDING |
| frozen model deployment | not yet authorized from 6-sample pilot alone | PENDING |

## Decision

Cohort F **passes the technical exact-representation gate** and is eligible to proceed to a frozen-model deployment pilot. It is not yet a primary external-validation cohort: subject-level independence, baseline/timepoint selection, and the effect of the fragment-length policy must be resolved before model scores are reported as validation evidence. No feature selection, threshold tuning, or model retraining was performed on Cohort F.
