# FBBP-DB Public Reproducibility Release

Version: `v20261119`

This repository provides the versioned, machine-readable reproducibility record for FBBP-DB.

Repository: https://github.com/changyufei222/fbbp-db-public-repro

Frozen version DOI: https://doi.org/10.5281/zenodo.24361542

Concept DOI for the latest version in this release line:
https://doi.org/10.5281/zenodo.20626952

## Included

- `Supplementary_Data_1_Field_Level_Data_Dictionary/`: field-level data dictionary.
- `Supplementary_Data_2_LLM_Prompts_and_Curation_QC/`: prompt templates, run metadata, G240 benchmark objects, and curation-QC metrics.
- `Supplementary_Data_3_Method_and_Release_Metadata/`: retrieval, structural, counting, and release metadata.
- `Supplementary_Data_4_Figure_Source_Data_Code_Manifest/`: figure source data and source-to-analysis mappings, including the reproducible PCSK9/4OV6 case-study data.
- `Supplementary_Data_5_Key_Source_Data/`: key statistical source tables and summaries.
- `Repository_Files_FBBP_DB_v1_release/`: database tables, schemas, code, manifests, test data, release documentation, and the complete frozen ML release.
- `MANIFEST_SHA256.csv`: checksums for every included file except the manifest itself.

## Scope

The release contains machine-readable data, source data, schemas, analysis and recomputation code, manifests, and technical documentation.

## Reproducibility Notes

The PCSK9/4OV6 source folder contains the exact D/E/G PDB subset, atom-pair contact table, per-residue ASA/ΔASA table, loop summary, recomputation script, and a standalone Figure S10D plotting script. The numerical scripts do not contain the manuscript-reported 17/180/68.1/524.8 values as constants.

## Licensing

Code is licensed under MIT (`LICENSE`). Curated data and documentation are licensed under CC BY 4.0 (`LICENSE-DATA`). Third-party source records remain subject to their original terms.
