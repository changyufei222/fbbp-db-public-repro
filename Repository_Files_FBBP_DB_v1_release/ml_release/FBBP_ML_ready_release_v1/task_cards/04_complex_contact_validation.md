# Task card: biological-assembly complex contact validation

## Dataset

- Binder-target contact tier: 53 records with copied deposited biological assemblies.
- Recognition-window target-contact tier: 19 records with the strictest contact evidence.

## Input and evidence fields

- Inputs: `assembly_coordinate_relpath`, `best_binder_chain_id`, `best_target_chain_id`, `binder_sequence`, `target_sequence`
- Evidence: `minimum_binder_target_distance_angstrom`, `minimum_window_target_distance_angstrom`, mapped-window count, and chain-match coverage.

## Recommended use

Use the 53-record tier for descriptive interface/contact analysis. Use the 19-record tier only as a high-confidence external validation or case-study set for contact/interface methods.

## Boundary

Neither tier provides enough independent examples for primary complex-model training. The positive-contact records do not provide matched experimental non-contact negatives.
