# GSE243474 exact 300-bp representation pilot v1

Two representative fragment BED files were processed against the archived A 300-bp stable feature universe (`A_stable_bins_top10pctCV_300bp.sorted.bed`). Fragments were retained when `20 <= end-start <= 1000` bp; longer or malformed records were excluded and not silently counted.

| sample | group | A features | nonzero features | nonzero fraction | fragment overlaps |
|---|---|---:|---:|---:|---:|
| GSM7788592 | Colorectal | 882,903 | 812,667 | 0.920 | 8,432,844 |
| GSM7788666 | Healthy | 882,903 | 716,321 | 0.811 | 3,157,118 |

Both samples produce finite counts on the exact A 300-bp hg19 feature universe without liftover. This is an assembly/representation feasibility pass, not a model-performance result. The remaining four pilot samples should be processed with the same fixed fragment policy before declaring Cohort F eligible for frozen deployment.
