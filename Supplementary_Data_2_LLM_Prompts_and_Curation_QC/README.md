# Supplementary Data 2

This folder contains the prompt templates, run metadata, benchmarking objects,
and sanitized curation/QC metrics used to evaluate LLM-assisted extraction.

## Benchmark sequence placeholder convention

`G240_benchmark_data_dictionary.tsv` records the interpretation of sequence
fields in the formal benchmark. In the G240 gold and prediction files, a
terminal run of `X` marks placeholder positions for which no concrete residue
was retained in the benchmark object. `X` is not imputed or treated as a
resolved amino-acid identity. Exact-sequence scoring compares the released
strings literally, so these terminal placeholders do not claim complete
end-to-end sequence fidelity. No `X` positions occur in the final FBBP-DB
protein sequences. Thirteen of the 240 gold records contain such terminal runs.

## Model identifiers

The working extraction pipeline records the provider model identifier
`gemini-3.1-pro-preview`. The formal G240 comparison is labelled
`Gemini 3.1 Pro`; its exact dated provider snapshot identifier was not retained
in the frozen run metadata. `DeepSeek-V3.2-Instruct` was used for the reviewer
and arbiter roles, and `Qwen3.5-122B-A10B` was used only in the comparison
configuration. These distinctions are preserved in
`Supplementary_Data_2_LLM_RUN_METADATA.tsv`.

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
- `G240_benchmark_data_dictionary.tsv` defines sequence placeholder handling
  and the exact-match interpretation for the formal benchmark.

## Prompt and QC files

The complete literature, patent, extractor, reviewer, and arbiter prompts are
provided under `prompt_files/`. Model-role assignments and benchmark provenance
are recorded in `Supplementary_Data_2_LLM_RUN_METADATA.tsv`. Sanitized
extractor-reviewer reliability metrics are under
`llm_curation_qc_metrics_sanitized/`.

The frozen role prompts are retained verbatim for provenance. Their internal
label `FBSP/FBSBP` is a legacy working acronym for the resource now released as
FBBP-DB; a terminology note at the start of each prompt makes this mapping
explicit.
