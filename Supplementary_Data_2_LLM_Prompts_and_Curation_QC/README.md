# Supplementary Data 2

This folder contains the prompt templates, run metadata, benchmarking objects,
and sanitized curation/QC metrics used to evaluate LLM-assisted extraction.

## Benchmark files

- `llm_curation_qc_metrics_sanitized/gold_benchmark_v1_runtime_ready.csv` is
  the formal 240-record G240 benchmark table used for source-level recovery
  and field-level evaluation.
  It contains 240 rows and the core scaffold, target, affinity, sequence,
  source, manual-verdict, and confidence fields used in the benchmark. The
  `Manual_Gold_Status` field contains the finalized review status (233 complete
  or corrected records and 7 partial-support records); the earlier preparation
  state is retained as `PreReview_Gold_Status`.
- `G240_builder_combined_arbiter_predictions.csv` contains the 289 combined
  prediction rows produced by the G240 multi-role workflow.
- `G240_benchmark_aligned_predictions.csv` contains the prediction rows aligned
  to the 240 gold records and the fields needed for field-level comparison.
- `G240_validity_row_report.csv` contains the record-level benchmark match
  summary.

## Prompt and QC files

The complete literature, patent, extractor, reviewer, and arbiter prompts are
provided under `prompt_files/`. Model-role assignments and benchmark provenance
are recorded in `Supplementary_Data_2_LLM_RUN_METADATA.tsv`. Sanitized
extractor-reviewer reliability metrics are under
`llm_curation_qc_metrics_sanitized/`.
