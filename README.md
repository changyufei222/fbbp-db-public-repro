# FBBP-DB Public Reproducibility Release

Version: `v20261119`

This package is the public, machine-readable reproducibility layer for FBBP-DB. It is deliberately separate from the private manuscript submission package.

## Included

- `Supplementary_Data_1_Field_Level_Data_Dictionary/`: field-level data dictionary.
- `Supplementary_Data_2_LLM_Prompts_and_Curation_QC/`: prompt templates, run metadata, G240 benchmark objects, and curation-QC metrics.
- `Supplementary_Data_3_Method_and_Release_Metadata/`: retrieval, structural, counting, and release metadata.
- `Supplementary_Data_4_Figure_Source_Data_Code_Manifest/`: figure source data and source-to-analysis mappings, including the reproducible PCSK9/4OV6 case-study data.
- `Supplementary_Data_5_Key_Source_Data/`: key statistical source tables and summaries.
- `Repository_Files_FBBP_DB_v1_release/`: database tables, schemas, code, manifests, test data, and release documentation.
- `MANIFEST_SHA256.csv`: checksums for every included file except the manifest itself.

## Deliberately Excluded

- Cover letter.
- Unsubmitted manuscript DOCX.
- Supplementary Information DOCX.
- Submission workbook containing publication-formatted supplementary tables.
- Rendered manuscript and supplementary figures.

These files remain in the private submission package and must not be added to a public release before the authors approve public dissemination.

## Reproducibility Notes

The PCSK9/4OV6 source folder contains the exact D/E/G PDB subset, atom-pair contact table, per-residue ASA/ΔASA table, loop summary, recomputation script, and a standalone Figure S10D plotting script. The numerical scripts do not contain the manuscript-reported 17/180/68.1/524.8 values as constants.

## Licensing

Code is licensed under MIT (`LICENSE`). Curated data and documentation are licensed under CC BY 4.0 (`LICENSE-DATA`). Third-party source records remain subject to their original terms.
