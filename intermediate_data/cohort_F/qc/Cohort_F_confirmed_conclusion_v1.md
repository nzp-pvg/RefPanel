# Cohort F confirmed conclusion v1

Cohort F passes technical exact frozen-representation transfer and frozen C-model deployment pilot gates. Six sample-level fragment BED files (3 Colorectal, 3 Healthy) were mapped to the frozen 20,000-feature hg19 model space without feature selection, threshold tuning, or retraining. The primary 20–1000 bp policy produced complete finite inputs and separated all 3 CRC samples above all 3 Healthy samples (pilot AUC 1.00; descriptive only). The 20–500 bp sensitivity policy gave essentially unchanged scores (CRC 0.934, 0.987, 0.990; Healthy 0.798, 0.760, 0.807), indicating the pilot signal is not driven by the excluded 500–1000 bp tail.

This is a technical/deployment pilot, not a definitive external-validation performance estimate. Formal cohort inclusion still requires subject-level clinical metadata and baseline/timepoint confirmation for the selected records.
