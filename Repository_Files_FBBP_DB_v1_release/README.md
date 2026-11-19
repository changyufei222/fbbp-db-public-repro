# FBBP-DB v20261119 data release

This directory is the reproducibility data package for the Briefings in Bioinformatics manuscript:

**Large language model-assisted Flexible-Backbone Binding Protein Database (FBBP-DB) and Analysis of Functional Recognition Regions**

Release date: 2026-11-19

## What is included

- `tables/`: normalized FBBP-DB schema/result tables, including protein-level records, structure context, active-loop annotations, flexibility evidence, predictions, audit, and provenance tables.
- `schema/`: SQL and Mermaid schema/ER assets.
- `sequence_structure/`: web-ready sequence/structure manifests and the combined FASTA for final protein records.
- `analysis/strict_background/`: strict AIR/ALL background comparison outputs and report.
- `analysis/air_all_background/`: enhanced loop-level AIR/ALL background comparison tables, summaries, and audit notes.
- `analysis/binding_interface/`: strict binding-window/interface residue tables, binding/non-binding statistics, the PCSK9/1459D05 case-study summary, and the complete 4OV6 contact/ASA recomputation folder.
- `figure_source_data/`: package-relative figure source-data/code manifests, schema/web-resource figure assets, and reproducible Figure S10D source-data plot generation.
- `code/`: code snippets/scripts located for strict-background, AIR/ALL, binding-interface and website-bundle analyses.
- `ml_release/`: the complete frozen binder-only and target-aware split objects, derived split views, structure-core coordinates, task cards, audit reports, and source scripts.
- `build_audit/`: sanitized public build-provenance summary without local absolute paths.
- `third_party_boundaries/`: redistribution policy for source records, public identifiers, and restricted full-text/raw third-party materials.

## Key current counts

- Protein-level records: `1996`
- Active-loop / recognition-window records: `3383`
- Loop-level flexibility evidence rows: `3383`

## Files to inspect first

1. `MANIFEST.tsv`
2. `DATA_DICTIONARY.tsv`
3. `schema/`
4. `tables/proteins.csv`
5. `tables/loop_annotations.csv`
6. `tables/loop_flexibility_public_summary.csv`
7. `analysis/ANALYSIS_ASSET_STATUS.md`
8. `build_audit/SANITIZED_BUILD_PROVENANCE.tsv`
9. `ml_release/ml_release_verified.json`

## Data dictionary

`DATA_DICTIONARY.tsv` is synchronized with the current field-level dictionary in `Supplementary_Data_1_Field_Level_Data_Dictionary/`.

## DOI and citation status

The version DOI for this frozen release is `10.5281/zenodo.24361542`. The concept DOI `10.5281/zenodo.20626952` points to the latest release in this line.

## Filename hygiene

Public package paths use English filenames. `filename_alias_map.csv` records the original non-ASCII working filenames and their English public-package paths.

## Public-path hygiene

Generated release documentation uses package-relative paths. `build_audit/SANITIZED_BUILD_PROVENANCE.tsv` provides release-path and header-preview provenance without local absolute source paths.
