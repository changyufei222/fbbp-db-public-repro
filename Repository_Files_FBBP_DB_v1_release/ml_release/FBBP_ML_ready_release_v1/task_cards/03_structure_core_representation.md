# Task card: structure-core representation and family prediction

## Dataset

- Mapped-95 tier: 201 experimental X-ray/cryo-EM records with copied PDB coordinates.
- Exact-sequence tier: 189-record subset with 100% retained-core sequence identity.

## Input and label

- Inputs: `coordinate_relpath`, `resolved_chain_id`, `binder_sequence`
- Available supervised label: `scaffold_category`
- Split: inherited frozen parent `fold_id` and `split_role`

## Recommended evaluation

Use the exact tier for strict structural benchmarking and the mapped-95 tier for broader robustness analysis. Report macro-F1 for family prediction or representation-transfer metrics, with the coordinate tier named in every result.

## Boundary

These PDB files are self-contained coordinate inputs, but the dataset is small for de novo 3D model training. It is appropriate for feature analysis, fine-tuning, or external-model evaluation, not large-scale structural pretraining.
