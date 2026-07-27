DATA EXPORT -- Boolean_2 correlation pipeline
Generated: 2026-07-28 00:40:02
======================================================================

FOLDER CONTENTS
----------------------------------------------------------------------
correlations_equal/
    Pearson & Spearman correlation tables (network_T_s,
    intra/inter_module_T_s, raw_module_T_s, impurity, density, ...
    vs module/network P_s and PCA metrics), equal (1/l_max)
    weighted influence matrix, for l_max = 1, 2, 5, 7, 10.

correlations_factorial/
    Same correlation tables, but using factorial (1/l!) weighted
    influence matrix, for l_max = 5, 7, 10.

module_level_data/equal/ and module_level_data/factorial/
    Per-module, per-network raw data points (network_T_s,
    module_T_s, module_pheno_score, impurity, density, etc.)
    used to compute the module-level correlation rows above --
    one file per l_max per attractor-selection type
    (dominant / top30 / top60).

panel1_data/
    panel1_data.csv            -- T_network, T_X_intra, T_X->Y
                                   for l_max = 5, 7, 10 (equal
                                   weighting), externally supplied.
    panel1_data_lfactorial.csv -- same, factorial weighting.
    (l_max = 1, 2 panel1-equivalent tables are built in-memory by
    the lmax1_2 scripts and are NOT saved to disk separately --
    their values are already folded into the module_level_data
    and correlations_equal CSVs above.)

pca_metrics/
    pca_metrics_by_folder.csv -- pc1_var and n_pcs_90 per network,
    used as the PC1_var / NumPC_90 correlation targets.

diagnostics/
    raw_vs_intra_diagnostic_lmax*.csv, if you generated them --
    compares raw_module_T_s (from Single_X.topo) against
    T_X_intra (from Comb_X_Y.topo) per module/network.

raw_module_T_s/
    raw_module_ts_{equal,factorial}_lmax{1,2,5,7,10}.csv
    Module strength (T_s) computed from each Single_X.topo file
    ALONE (i.e. the module in isolation, not embedded in a
    combined network), per module, per l_max, per weighting.
    Previously only lived in memory during a correlation run.

network_level_P_s/
    network_ps_{equal,factorial}_lmax{1,2,5,7,10}.csv
    Whole-network phenotypic score (not per-module), per network,
    per state_type (dominant/top30/top60), per l_max, per
    weighting. Includes a p_s_formula column with BOTH formulas
    computed side by side:
      'state_based'        -- signed: team1_mean_state - team2_mean_state
      'influence_weighted' -- unsigned: abs(team1_input - team2_input)
    Previously only lived in memory during a correlation run.

panel1_equivalent_lmax1_2/
    panel1_equivalent_{equal,factorial}_lmax{1,2}.csv
    T_network, T_X_intra, T_X->Y per network, for l_max = 1, 2
    only (l_max = 5,7,10 already covered by panel1_data.csv /
    panel1_data_lfactorial.csv above). The factorial l_max=1,2
    version did not exist in any prior script -- newly computed
    here using the same factorial-weighted influence matrix as
    the l_max=5,7,10 factorial pipeline.

NOTES ON PHENOTYPIC SCORE (P_s) VERSIONS
----------------------------------------------------------------------
Two different P_s formulas were used at different points; check
which script/date produced a given correlations_*.csv if it
matters for your analysis:
  - State-based (signed):   team1_mean_state - team2_mean_state
  - Influence-weighted (unsigned): abs(team1_input - team2_input),
    where team_input = mean_i( sum_j Inf[i,j] * state[j] )

LOG OF THIS EXPORT
----------------------------------------------------------------------
  OK    D:\csb_lab\Boolean_2\correlations_pearson_equal_lmax1.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_equal_lmax2.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_equal_lmax5.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_equal_lmax7.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_equal_lmax10.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_equal_lmax1.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_equal_lmax2.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_equal_lmax5.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_equal_lmax7.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_equal_lmax10.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_factorial_lmax5.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_factorial_lmax7.csv
  OK    D:\csb_lab\Boolean_2\correlations_pearson_factorial_lmax10.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_factorial_lmax5.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_factorial_lmax7.csv
  OK    D:\csb_lab\Boolean_2\correlations_spearman_factorial_lmax10.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax10_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax10_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax10_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax1_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax1_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax1_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax2_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax2_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax2_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax5_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax5_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax5_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax7_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax7_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_equal\lmax7_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax10_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax10_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax10_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax5_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax5_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax5_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax7_dominant_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax7_top30_module_data.csv
  OK    D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\lmax7_top60_module_data.csv
  OK    D:\csb_lab\Boolean_2\panel1_data.csv
  OK    D:\csb_lab\Boolean_2\panel1_data_lfactorial.csv
  OK    D:\csb_lab\Boolean_2\pca_plots\pca_metrics_by_folder.csv
  MISS  D:\csb_lab\Boolean_2\Team_strength_pearson_equal\raw_vs_intra_diagnostic_lmax*.csv  (no files matched -- skipped)
  MISS  D:\csb_lab\Boolean_2\Team_strength_pearson_factorial\raw_vs_intra_diagnostic_lmax*.csv  (no files matched -- skipped)
  MISS  D:\csb_lab\Boolean_2\Raw_data_supplement\raw_module_T_s\raw_module_ts_*.csv  (no files matched -- skipped)
  MISS  D:\csb_lab\Boolean_2\Raw_data_supplement\network_level_P_s\network_ps_*.csv  (no files matched -- skipped)
  MISS  D:\csb_lab\Boolean_2\Raw_data_supplement\panel1_equivalent_lmax1_2\panel1_equivalent_*.csv  (no files matched -- skipped)
