# Task card: scaffold-family classification

## Dataset

`binder_full_length_ml_ready_8_1_1_v1.csv` contains 1,996 full-length binder sequences across 12 normalized scaffold families. The frozen 8:1:1 split contains 1,595 train, 199 validation, and 202 test objects.

## Input and label

- Input: `binder_sequence`
- Label: `supervised_label` / `scaffold_category`
- Leakage control: 95% within-family global-affine binder clusters never cross fold boundaries.

## Recommended evaluation

Report accuracy, macro-F1, per-family recall, confusion matrix, and cluster-aware split provenance. Macro-F1 is required because family sizes are uneven.

## Boundary

This is a scaffold-family classification task. It is not a binding/non-binding prediction task and does not establish affinity prediction.
