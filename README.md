# 🧬 Frozen cross-cohort plasma cfMeDIP-seq transfer hub

This repository provides a reproducible analysis workflow for testing whether a plasma cfMeDIP-seq representation can be transferred across cohorts without target-cohort feature selection or model refitting. The framework follows a transparent **Tune → Calibrate → Freeze → Deploy** design. Cohort A defines a stable reference representation. Cohort B calibrates input-control reliability rules. The resulting feature order, preprocessing configuration and model are frozen before deployment to independent cohorts C, D and F.

The repository is organized for traceability, parameter transparency and figure-level reruns. It separates frozen artifacts from external source data so that a reviewer can inspect the locked representation and deployment results without downloading the complete fragment archive.

---

## 🎯 What you can reproduce here

- ✅ Corrected Cohort A stability analysis across 300, 600, 900 and 1200 bp windows
- ✅ Cohort B paired IP/input-control concordance and IC-quantile sensitivity analysis
- ✅ A/B-derived frozen 300-bp feature order with 20,000 model features
- ✅ Frozen model deployment in Cohort C SCLC data
- ✅ Cohort D HNSCC cfDNA primary deployment with paired PBL secondary analysis
- ✅ Conventional, D-adaptive and A/B-locked strategy benchmarking in Cohort D
- ✅ Cohort D cross-material correlation and leukocyte-background QC
- ✅ Cohort F cross-center multi-cancer frozen-representation deployment summaries
- ✅ Explicit deployment failure flags and metadata-integrity audits
- ✅ PNG figure-generation scripts for the integrated C/D/F deployment figure

---

## 🚀 Quick start

### 1) Environment

The analysis uses R with `data.table`, `glmnet`, `matrixStats`, `dplyr`, `stringr`, `readr`, `tidyr`, `ggplot2` and `patchwork` as applicable. Python audit scripts use the standard library. Fragment-to-feature construction requires `bedtools`.

### 2) Recommended entry points

Run the compact deployment checks from the repository root:

```bash
Rscript code/analysis/deploy_frozen_model_C_D_v1.R
Rscript code/analysis/benchmark_D_strategies_v1.R
Rscript code/analysis/build_transfer_benchmark_CDF_v1.R
Rscript code/figures/build_Figure6_expanded_v2.R
```

The first command regenerates Cohort D scores from the included model-space matrices. Cohort C can be regenerated only when its authorized source matrices are supplied with `--c-cfdna`, `--c-ncc` and `--c-meta`.

For Cohort F raw fragment processing, download the BED files listed in `intermediate_data/cohort_F/source/GSE243474_MeDIP_bed_download_manifest_v1.tsv`, then run:

```bash
bash code/analysis/build_cohort_F_frozen_representation_v1.sh /path/to/GSE243474_bed_directory
Rscript code/analysis/score_cohort_F_frozen_model_all_v1.R
Rscript code/analysis/bootstrap_cohort_F_deployment_v1.R
```

The public release does not contain the 35 GB BED archive. The script stops with an explicit input message when the directory is absent.

---

## 🧭 Study design and cohort roles

| Cohort | Role | Samples or records | Primary use | Boundary |
|---|---|---:|---|---|
| **A** | Healthy reference and development panel | 35 healthy controls; 90 HCC interpretation samples | Stable genomic representation | No target-cohort feature selection |
| **B** | Transfer calibration panel | 25 matched IP/IC pairs | Input-control gate and representation reliability | Calibration only |
| **C** | Independent SCLC deployment | 94 cfDNA samples; 74 paired PBL samples | Locked cancer deployment | No C-specific feature selection in the locked analysis |
| **D** | Independent HNSCC stress test | 50 cfDNA; 49 eligible PBL | Primary cfDNA transfer and secondary PBL background analysis | `PBL_control_20` excluded after byte-identical FaDu audit |
| **F** | Cross-center multi-cancer portability | 301 MeDIP records; 283 retained title-derived baseline proxy units | Frozen representation feasibility and cancer-stratified deployment | Subject-level metadata remains incomplete for definitive clinical validation |
| **E** | Screened but excluded | 98 promoter-space records | Integrity and compatibility audit | Promoter-space mismatch and unresolved mixed biological units |

### Cohort D integrity rule

The primary D cfDNA analysis retains 30 HNSCC samples and 20 healthy controls. The paired PBL analysis retains 30 HNSCC samples and 19 healthy controls. `PBL_control_20.bed` was byte-identical to the deposited FaDu profile and is excluded only from PBL analyses. The corresponding `cfDNA_control_20` file remains eligible for cfDNA-only analysis.

### Cohort F interpretation rule

Cohort F contains 34 healthy proxy units and 249 cancer proxy units across 15 cancer groups after the current title-based deduplication. These records demonstrate technical portability in the frozen feature space. They are not equivalent to a fully adjudicated independent-subject clinical cohort until subject identity, baseline status and repeated sampling are resolved from complete metadata.

---

## 🔒 Frozen protocol

### Reference representation

- Primary genomic resolution: **300 bp**
- Stability candidate rule: lowest **10% raw-count CV** in the healthy reference panel
- Blacklist handling: whole windows are excluded when any component overlaps a blacklisted interval
- Multi-scale sensitivity: 300, 600, 900 and 1200 bp windows are constructed from the complete grid before blacklist filtering
- Low-signal protection: low-mean or zero-signal windows are excluded from stability selection

### Input-control calibration

- Primary IC gate: **20th percentile rule** from Cohort B
- Archived MAIN threshold: `1.98328` on the B mean-centered IC scale
- Ratio pseudocount: **0.1**
- IC data are used for calibration only. Cohorts C, D and F are deployed without requiring new input-control libraries.

### Frozen model artifact

The authoritative files are in `intermediate_data/frozen_model/`:

- `C_model_features_v1.tsv`: 20,000 ordered feature coordinates
- `C_model_features_v1.bed`: corresponding interval table
- `C_frozen_model_v1.rds`: serialized fitted model and feature order
- `C_frozen_model_v1_manifest.tsv`: model provenance and fixed transformation

The archived model uses a binomial elastic-net model with `alpha = 1`, standardization enabled, the recorded `lambda.1se` and the fixed transformation `log2(count_or_signal + 1)`. The model was trained on 124 A samples comprising 89 HCC samples and 35 controls. Cohort D and Cohort F are not used to alter this feature order or these coefficients.

---

## 📊 Evaluation outputs

The compact deployment tables are stored in `intermediate_data/deployment/`.

- Cohort C frozen cfDNA deployment: AUC **0.864** in the current released summary
- Cohort D frozen cfDNA deployment: AUC **0.663** for 30 HNSCC versus 20 healthy controls
- Cohort D frozen PBL secondary deployment: AUC **0.507** for 30 HNSCC versus 19 eligible healthy controls
- Cohort D strategy benchmark: conventional, D-label-adaptive and A/B-locked scores are kept together. The adaptive result is a benchmark only and is not a replacement for the locked endpoint.
- Cohort F: cancer-group AUC estimates with stratified bootstrap intervals are retained as portability evidence. Near-chance groups remain visible as deployment boundaries.

The D PBL result is therefore reported as a secondary failure-boundary signal. It is not silently combined with the primary cfDNA endpoint.

---

## 🧪 Quality-control and leakage safeguards

- No C, D or F labels are used to choose the frozen feature order.
- No D or F model is refit for the locked endpoint.
- Paired PBL analysis is labeled separately from cfDNA deployment.
- `PBL_control_20` is retained as an integrity anomaly record even though it is excluded from eligible PBL scoring.
- Cohort F title-proxy deduplication is recorded in `GSE243474_subject_audit_v1.tsv` and `GSE243474_locked_manifest_v1.tsv`.
- Explicit failure flags are written to `D_deployment_failure_flags_v1.tsv`.
- The audit in `intermediate_data/provenance/ABCD_pipeline_audit_v1.md` documents legacy leakage risks and the corrected lock boundary.

---

## 📁 Repository layout

```text
MS_5_code_id/
├── README.md
├── code/
│   ├── analysis/
│   └── figures/
├── intermediate_data/
│   ├── frozen_model/
│   ├── cohort_A/
│   ├── cohort_B/
│   ├── cohort_C/
│   ├── cohort_D/derived_v1/
│   ├── cohort_F/{derived,qc,source}/
│   ├── deployment/
│   └── provenance/
└── .gitignore
```

---

## 🌍 External source data

The release does not redistribute source archives. Use the original records for raw-data reconstruction:

- Cohort A: Zenodo **11251606**
- Cohort B: GEO **GSE152631**
- Cohort C: Zenodo **7235989**
- Cohort D: Zenodo **4698004**
- Cohort F: GEO **GSE243474**

The compact release is sufficient for frozen-model inspection, D deployment regeneration, QC review and figure reconstruction. Complete raw-to-result regeneration additionally requires the source matrices or fragment files excluded above.

---

## ♻️ Reproducibility and reporting guardrails

This repository separates three evidence levels:

1. **Frozen deployment:** a model and feature order fixed before the target cohort is scored.
2. **Technical portability:** an external sample can be mapped into the frozen representation with finite values and no feature reselection.
3. **Clinical validation:** an independently adjudicated subject-level cohort with complete baseline and outcome metadata.

Cohort F currently supports the second level. It should not be described as definitive clinical validation until the metadata gap is resolved. Similarly, the D PBL result is a stress-test boundary for leukocyte-background portability rather than a primary cancer-classification endpoint.

---

## 📜 Citation

When using this repository, cite the accompanying manuscript and the original data records listed above. Software versions, fixed parameters, model artifacts and intermediate result tables should be reported together so that the locked deployment can be distinguished from exploratory adaptive benchmarks.
