# ML-ready data dictionary

## Common fields

| Field | Meaning |
|---|---|
| `resource_record_id` | Stable FBBP parent protein identifier |
| `sequence_object_id` | Stable normalized binder-sequence identifier |
| `scaffold_category` | Normalized FBBP family label |
| `binder_sequence` | Full-length released binder amino-acid sequence |
| `fold_id` | Frozen 10-fold assignment unit |
| `split_role` | `train`, `validation`, or `test` in the fixed 8:1:1 view |

## Binder table

`binder_full_length_ml_ready_8_1_1_v1.csv` has 1,996 rows. `supervised_label` is the normalized scaffold family, so it supports a 12-class scaffold-family classification task. It does not provide binding/non-binding labels.

## Target-pair table

`target_sequence_aware_pair_ml_ready_8_1_1_v1.csv` has 710 unique protein-target pairs from 679 unique binder objects and 145 frozen target sequences. `curated_target_association_label=1` means that the pair is a curated sequence-resolved association in FBBP-DB. It is not a verified negative-pair dataset: all unobserved pairs remain unknown, not negative.

`fasta_file` values are relative to the `ml_release` package root and resolve under `FBBP_ML_target_sequence_aware_10fold_v1/data/fasta/`.

## Structure-core tables

`structure_core_mapped_95_ml_ready_8_1_1_v1.csv` has 201 records, and `structure_core_exact_sequence_ml_ready_8_1_1_v1.csv` has the 189-record exact subset. `coordinate_relpath` is relative to the package root and points to a copied PDB file. `coordinate_sha256`, `parsed_atom_count`, and `parsed_chain_ids` verify the coordinate payload.

## Complex tables

The 53-record binder-contact and 19-record recognition-window-contact tables contain `assembly_coordinate_relpath`, pointing to gzip-compressed PDB biological assemblies. They include sequence-resolved binder and target chain IDs, target accession, contact distance, and inherited 8:1:1 split roles. The 19-record tier is a high-confidence validation tier, not a training set.
