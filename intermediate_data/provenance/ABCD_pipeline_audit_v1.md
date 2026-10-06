# A→B→C/D frozen-transfer audit (v1)

## Evidence-based reconstruction

| Stage | What the archived code actually does | Audit conclusion |
|---|---|---|
| A | Uses 35 CTL samples after excluding N19; removes blacklist IDs; creates 300/600/900/1200-bp matrices; calls the lowest 10% raw-count CV bins stable. | A defines healthy-control low-variability candidate reference regions, not classifier features. Existing larger bins were formed after blacklist deletion and can bridge genomic gaps. Low/zero-signal bins can also appear artificially stable. |
| B | Uses 25 paired A/C keys with M and IC tracks; compares `log2(M+1)` with `log2((M+pc)/(IC+pc))`; sweeps IC median gates and records q=0.20, pc=0.1 as MAIN. | B calibrates an input-control preprocessing setting on the 300-bp blacklist-filtered universe. It does not select a genomic scale from A and does not emit locked feature IDs. The archived MAIN QC contains 8,829,028 bins before and 7,063,220 after the q20 IC gate. |
| C | Intersects C matrices with the A 300-bp blacklist-filtered BED, applies `log2(count+1)`, computes paired SCLC cfDNA−PBL, but subtracts a pooled SCLC PBL reference from NCC. Later blocks select the 20,000 highest-variance bins from C before PCA/CV. | C contains target-cohort feature selection and therefore is not a frozen external transfer analysis. PBL preprocessing is asymmetric between SCLC and NCC. The monolithic script also has hard-coded CVD_MS_6 paths and depends on in-memory objects created by earlier blocks. |
| D | Contains 30 HNSCC + 20 healthy cfDNA and 30 HNSCC + 20 healthy PBL BED files. PBL_control_20 is byte-identical to FaDu_MeDIP. | D must remain feature-selection-free. Primary validation is 30+20 cfDNA. Secondary PBL analysis is 30+19 after excluding PBL_control_20; cfDNA_control_20 remains valid for cfDNA-only analysis. |

## Frozen design implemented

A v2 full 300-bp grid → construct truly contiguous scales → blacklist whole windows if any component is blacklisted → exclude bins below an explicit mean-signal floor → define CTL low-CV candidates → quantify sample-subsampling stability → express cross-scale overlap on a common 300-bp atomic coverage universe. The primary lock is then fixed at 300 bp and uses A stability plus the B-only IC q20 detectability gate; C and D are explicitly excluded from feature selection. C/D may only be subset and ordered by `AB_locked_features_v2.tsv`.

The A and D interval tables use 1-based inclusive coordinates (`chr1 1 300`, `chr1 301 600`) even though filenames use `.bed`. They must be documented and parsed as 1-based inclusive interval tables, not silently converted under standard 0-based half-open BED assumptions. Cohort D raw filenames and 27-GB source files remain immutable provenance; derived manifests, QC and compact locked matrices are stored separately.

## Leakage and preprocessing boundary

No variance ranking, differential testing, lasso screening, PCA loading selection, threshold tuning, or missing-feature-dependent reselection may be performed in C or D. A model intended for direct C/D comparison must use the same analyte transformation for cases and controls. The prior C paired-patient-PBL versus pooled-patient-PBL-for-controls contrast is retained only as exploratory asymmetric preprocessing, not primary validation. For D, cfDNA-only is the primary symmetric analysis; paired PBL is a separately labelled secondary analysis and cannot silently replace the primary endpoint.

## Reproducibility status and stop condition

The repository snapshot does not contain the large A count/sample inputs, A stable-ID files, B M/IC matrices or executable B avg files, nor the C count matrices. Therefore the corrected scripts and D manifest can be produced and syntax/QC checked now, but a truthful final locked feature artifact and compact C/D matrices cannot be materialized until A/B inputs or a previously generated locked-ID artifact are restored. Summary TSVs are insufficient to reconstruct genomic feature identities. This is an exact data-availability stop, not a reason to derive features from C or D.
