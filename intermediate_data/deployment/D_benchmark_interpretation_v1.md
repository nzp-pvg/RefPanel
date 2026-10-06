# Cohort D strategy benchmark and deployment flags (v1)

All strategies use the same 20,000-feature D matrices. `locked_A_model` is the only external-transfer result: it is trained in A and scored in D without D labels. `conventional` is the pre-specified mean log2 signal summary. `adaptive_5fold` is a D-internal comparator and is not external validation.

| Material | Conventional AUC | Adaptive 5-fold AUC | Locked A-model AUC | Decision |
|---|---:|---:|---:|---|
| cfDNA | 0.4933 | 0.4000 | 0.6633 | Retain as weak exploratory portability signal; no improvement claim. |
| PBL | 0.6105 | 0.3772 | 0.5070 | Deployment failure flag; do not report as reliable discrimination. |

The locked model shows some cfDNA portability above the simple global-signal comparator, but the absolute AUC is modest and the cohort is small and cross-disease. The PBL result fails the pre-specified discrimination gate (<0.60) and should be interpreted only as a negative secondary background-adjustment test. Adaptive scores are unstable in this small cohort and use D labels; they are shown only to demonstrate the cost of cohort-specific adaptation, not as a replacement for the locked result.

No feature retuning, threshold tuning, or label-informed harmonization was used for the locked results. The current evidence supports a methods-framework/stress-test claim, not a clinical classifier claim or a demonstrated performance improvement.
