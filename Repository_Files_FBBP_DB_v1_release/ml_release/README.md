# FBBP-DB frozen ML release

This directory contains the complete frozen machine-learning release described
in the Briefings in Bioinformatics manuscript.

## Contents

- `FBBP_ML_binder_full_length_10fold_v1/`: binder-only sequence clusters,
  ten-fold assignments, family-balance reports, and leakage audits.
- `FBBP_ML_target_sequence_aware_10fold_v1/`: frozen target sequences, target
  clusters, connected components, ten-fold assignments, and audit reports.
- `FBBP_ML_ready_release_v1/`: ML-ready binder, target-pair, structure-core,
  and complex-contact objects with task cards.
- `FBBP_ML_derived_split_views_20260716/`: deterministic 9:1, 8:2, 8:1:1, and
  7:3 views derived from the frozen ten-fold assignments.
- `source_scripts/`: scripts used to construct and audit the released splits.

## Integrity

`MANIFEST_SHA256.csv` lists every file in this directory except the manifest
itself. `ml_release_verified.json` records the primary object counts and split
crossing checks. Verification can be repeated from the public repository root:

```bash
python Repository_Files_FBBP_DB_v1_release/scripts/validate_release.py
```
