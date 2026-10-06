# GSE243474 BED pilot QC v1

Six sample-level fragment BED files were downloaded: three Healthy and three Colorectal MeDIP samples. All files are six-column BED records with canonical `chr` chromosome names. The pilot files are therefore coordinate-compatible in assembly naming with the existing A/B hg19 grid at the basic syntax level.

| group | files | approximate records | records >1 kb | interpretation |
|---|---:|---:|---:|---|
| Healthy | 3 | 4.3M each (range) | ~0.06% in inspected file | mostly short fragments |
| Colorectal | 3 | 13.8–15.5M | ~0.45–0.57% | mostly short fragments but non-negligible long-span tail |

The BED records include a small tail with very large spans (including >1 Mb). These records must not be silently treated as ordinary cfMeDIP fragments. The exact representation builder must define and report a fragment-length policy (for example, retain valid short fragments and separately quantify excluded long-span records), then compare sensitivity to that policy. This pilot is not yet an exact frozen-transfer pass; it is a successful assembly/format feasibility check with a pending fragment-integrity decision.
