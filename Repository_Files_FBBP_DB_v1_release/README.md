# FBBP-DB v20261119 data release

This directory is the reproducibility data package for the Briefings in Bioinformatics manuscript:

**FBBP-DB: a family-aware active-loop annotation framework for comparable recognition-region analysis of non-antibody scaffold binders**

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

## Data dictionary

`DATA_DICTIONARY.tsv` is synchronized with the current field-level dictionary in `Supplementary_Data_1_Field_Level_Data_Dictionary/`.

## DOI and citation status

The version-resolving concept DOI for this release line is `10.5281/zenodo.20626952`.

## Filename hygiene

Public package paths use English filenames. `filename_alias_map.csv` records the original non-ASCII working filenames and their English public-package paths.

## Public-path hygiene

Generated release documentation uses package-relative paths. `build_audit/SANITIZED_BUILD_PROVENANCE.tsv` provides release-path and header-preview provenance without local absolute source paths.
