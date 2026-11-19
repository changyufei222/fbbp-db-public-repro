# Task card: target-association retrieval

## Dataset

`target_sequence_aware_pair_ml_ready_8_1_1_v1.csv` contains 710 sequence-resolved curated protein-target pairs from 679 unique binder objects and 145 frozen target sequences. The split is controlled jointly by 95% binder clusters and 40% target-sequence clusters.

## Input and label

- Inputs: `binder_sequence`, `target_sequence`
- Positive association indicator: `curated_target_association_label=1`
- Split fields: `fold_id`, `split_role`, `target_cluster_40_id`, and `target_aware_component_id`

## Recommended evaluation

Use retrieval or ranking metrics such as Recall@K, MRR, and median rank against a predefined candidate target panel restricted to the same split role. Record the candidate-panel construction exactly.

## Boundary

No experimentally verified negative pairs are supplied. Unobserved pairs are unknown, not negatives. Do not train a binary binding classifier by treating random unobserved protein-target pairs as confirmed negatives. The target-aware layer is adequate for overall evaluation but dominated by adnectin and OBody; report family composition.
