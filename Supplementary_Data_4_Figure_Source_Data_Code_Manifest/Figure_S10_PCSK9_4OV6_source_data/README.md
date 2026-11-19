# PCSK9/4OV6 source data

`4OV6_PCSK9_Adnectin_D_E_G.pdb` is a local subset of PDB 4OV6 containing PCSK9 chains D/E and Adnectin chain G. The coordinate file is used only as the fixed structure input; water and other non-polymer records are ignored by the analysis.

Recompute the source tables from this directory:

```bash
python -m pip install -r requirements.txt
python recompute_pcsk9_1459D05_interface.py \
  --pdb 4OV6_PCSK9_Adnectin_D_E_G.pdb \
  --outdir .
```

The reported tables were generated with Biopython 1.88.

Definitions:

- Binder: Adnectin chain G.
- Target: PCSK9 chains D and E.
- Contacts: non-hydrogen atom pairs crossing the binder-target interface at a 4.5 Å cutoff. This is distinct from the unique residue-pair metric used for database-wide FRR contact summaries.
- ΔASA: per-residue solvent-accessible surface area in chain G alone minus chain G in the G-D/E complex, calculated with Biopython Shrake-Rupley using 960 surface points and a 1.4 Å probe radius.
- Only standard polymer residues are analyzed. Binder contact-residue counts distinguish insertion codes.
- `pcsk9_1459D05_loop_interface_summary.tsv` also reports the mean binder-alone ASA, bound-complex ASA, and ΔASA across all residues in each loop.

Reference loop-level values:

| FRR | Atom-pair contacts | ΔASA (Å²) |
|---|---:|---:|
| L2 (22-30) | 17 | 68.099 |
| L6 (60-66) | 0 | 0.000 |
| L7 (76-83) | 180 | 524.777 |

Source files:

- `pcsk9_1459D05_heavy_atom_contacts.tsv`: all 324 non-hydrogen atom-pair contacts.
- `pcsk9_1459D05_residue_burial_deltaASA.tsv`: per-residue contact and ASA/ΔASA values.
- `pcsk9_1459D05_loop_interface_summary.tsv`: aggregate values for L2, L6, L7, all FRRs, and all binder-interface residues.
- `pcsk9_1459D05_case_study_summary.tsv`: updated curated loop metadata and status fields.
- `recompute_pcsk9_1459D05_interface.py`: cleaning from raw structure to the three computed TSV files.
- `plot_Figure_S10D_pcsk9_deltaASA.py`: regenerates the PCSK9/4OV6 panel-D bar plot from the loop summary TSV without manuscript-value constants.
- `4OV6_PCSK9_Adnectin_D_E_G.pdb`: exact coordinate subset used for all calculations.
