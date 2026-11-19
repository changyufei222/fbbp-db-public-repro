# Figure source data for Supplementary Figures S2, S4, S5, and S6

This folder contains the loop-level and background tables needed to reproduce
the flexibility enrichment and avgRSA-ITsFlexible association analyses.

- `loop_flexibility_public_summary.csv`: the 3,383-FRR public summary table.
  The S4 association analysis uses the 3,301 FRRs with valid ITsFlexible
  scores and avgRSA values.
- `tbp_vs_airs_air_reference_v1.csv`: the 21,104-loop ALL-conformations
  antibody/TCR CDR background.
- `tbp_vs_airs_length_standardized_summary_v2.csv`: length-standardized
  observed-versus-expected fractions.
- `tbp_vs_airs_logit_summary_v2.csv`: length-adjusted logistic models.
- `tbp_vs_airs_mantel_haenszel_summary_v2.csv`: length-stratified sensitivity
  models.
- `Supplementary_Figure_S6_binding_leave_one_family_out.csv`: figure source data
  for the family-exclusion sensitivity analysis in Supplementary Figure S6.
- `strict_combined_background_summary_20260531.csv`: summary of the stricter
  combined background.
- `strict_combined_background_air_clusters_20260531.csv`,
  `strict_combined_background_fbbp_clusters_20260531.csv`, and
  `strict_combined_background_cluster_members_20260531.csv`: source tables for
  the strict-background comparison.

The protein-cluster bootstrap confidence interval shown in Figure S4 is
reported in
`../../Supplementary_Data_5_Key_Source_Data/Key_Result_Summaries/clustered_inference_summary.json`.
