# FBBP target-sequence-aware 10-fold audit

Status: **PASS - frozen target sequences and 40% target-aware fold IDs**

- release_name: target_sequence_aware_10fold_v1
- status: PASS
- release_date: 2026-07-16
- primary_target_identity_threshold: 0.4
- sensitivity_versions: [0.3, 0.5, 'exact_accession']
- frozen_target_sequences: 145
- target_sequence_resolved_protein_objects: 679
- target_aware_component_count: 123
- largest_target_aware_component: 42
- records_excluded_from_target_sequence_graph: 112
- name_target_records_excluded: 107
- accession_records_excluded_sequence_unavailable: 5
- uniparc_recovered_target_sequences: 4
- fold_sizes: {'0': 68, '1': 68, '2': 68, '3': 68, '4': 68, '5': 68, '6': 68, '7': 68, '8': 67, '9': 68}
- binder_cluster_crossings: 0
- target_cluster_40_crossings: 0
- component_crossings: 0
- derived_train_validation_test_view: FBBP_ML_derived_split_views_20260716
- Unresolved or ambiguous name-only target records are excluded rather than inferred.
- No train/validation/test CSV or benchmark ZIP is created.
