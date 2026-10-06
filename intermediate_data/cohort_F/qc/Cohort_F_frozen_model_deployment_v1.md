# Cohort F frozen C-model deployment assessment v1

The six pilot samples were converted into the A 300-bp stable feature universe and then checked against the frozen C model feature order. The frozen C model requires 20,000 model-specific feature IDs, whereas the pilot matrices currently contain the 882,903 A-stable-bin universe. Only **11/20,000** frozen model feature IDs are present in the current pilot count files.

Consequently, the preliminary score vector is not valid for deployment: zero-filling the missing 19,989 model features would create artificial near-constant predictions. No performance claim is made from those scores.

| gate | result |
|---|---|
| exact hg19 300-bp representation | PASS (6/6) |
| fixed fragment policy | PASS (20–1000 bp; start >= 0) |
| frozen C model input compatibility | **FAIL/PENDING** (11/20,000 features overlap) |
| independent subject/timepoint eligibility | PENDING |
| formal frozen deployment score | **NOT AUTHORIZED** |

The correct next action is to rebuild the six pilot matrices directly on the frozen C model's 20,000 feature BED, using the same fragment policy, then rerun prediction. This is a representation extraction issue, not evidence that Cohort F lacks biological portability.
