# Combined module-level dataset

`combined_module_level_data.csv` — 1440 rows x 29 columns, built by
`build_combined_dataset.py` (rerun it to regenerate).

## Grain

One row per **(pheno_score_method, l_max, network_name, module_abbr)**.

2 methods x 4 l_max x 56 networks, where the 56 networks are all module
combinations of the 6 modules AS, C, E, M, RS, PS:

| n_modules | networks | module rows per method/l_max |
|---|---|---|
| 2 | 15 = C(6,2) | 30 |
| 3 | 20 = C(6,3) | 60 |
| 4 | 15 = C(6,4) | 60 |
| 5 |  6 = C(6,5) | 30 |

## Columns

**Keys / identifiers**

| column | meaning |
|---|---|
| `pheno_score_method` | `newNorm` (the unsuffixed `correlations_*` folders) or `stateps` (the `_stateps` folders). A separate phenotypic-score calculation — keep them separate; do not pool. |
| `l_max` | 1, 5, 7, 10 |
| `n_modules` | 2-5 |
| `network_name` | e.g. `Comb_AS_C_E` |
| `modules_str` | module list in source order, from the panel file |
| `modules_sorted` | same list alphabetised, for joining/grouping |
| `has_PS` | 1 if the network contains the PS module. `has_PS == 0` exactly reproduces the `_noPS` folders (verified). |
| `module_abbr` | the module this row describes |
| `point_label` | original plot label, e.g. `AS(AS_C_E)` |

**Module-level**

| column | meaning |
|---|---|
| `module_T_s` | team strength of this module (identical to the panel's `T_<module>_intra`) |
| `P_s_dominant`, `P_s_top30`, `P_s_top60` | phenotypic score under the three state-selection variants. The three source files (`*_dominant/top30/top60_module_data.csv`) are identical except for this value, so they are pivoted into columns here. |
| `T_module_out_sum` | sum of directed strengths from this module to the network's other modules |
| `T_module_in_sum` | sum of directed strengths into this module |

**Network-level** (repeated across the network's module rows)

| column | meaning |
|---|---|
| `network_T_s` | network team strength (panel `T_network`) |
| `impurity`, `density`, `intra_density`, `inter_density` | network topology metrics |
| `T_network_intra_sum` | sum of all module `T_*_intra` (equals the sum of `module_T_s`) |
| `T_network_inter_sum` | sum of all directed inter-module strengths |

**PCA** (per network, l_max-independent; from `<Tag>_PCA_variance.csv`)

| column | meaning |
|---|---|
| `n_nodes`, `n_samples` | network size and RACIPE sample count |
| `n_PCs` | non-missing PCs in the spectrum |
| `PC1_var` | variance explained by PC1 |
| `NumPC_90` | PCs needed to reach 90% cumulative variance |
| `PCA_entropy` | Shannon entropy of the normalised variance spectrum |
| `PCA_entropy_norm` | `PCA_entropy / log(n_PCs)` |

## Validation

- Every network in the correlations folders has a matching `*_PCA_variance.csv`
  row, so the PCA join is complete (no missing values anywhere in the file).
- The four PCA metrics are **derived**, not copied. Their definitions were
  chosen to reproduce the published `correlations_pearson_equal_*.csv`
  coefficients: **1120/1120 correlation checks** across both methods, all
  l_max, all n_modules, all 7 metric rows and 5 network columns agree to
  <6e-5 (the rounding in the published files).
- `newNorm` vs `stateps` differ only in the `P_s_*` columns; the `T_s` columns
  differ by at most 4.4e-16 (summation-order noise).

## Known gaps

1. **Network-level P_s is not included.** The `dominant_network`, `top30_network`
   and `top60_network` rows of the correlations tables are network-level
   phenotypic scores, but only the *correlation coefficients* were saved — the
   per-network values are not stored anywhere for these networks, and they are
   not recoverable from the module scores (mean/max/min/sum all fail to
   reproduce the published coefficients).
2. **Directed pairwise strengths are summarised, not itemised.** The panel files
   hold a full directed matrix `T_<A>→<B>`; at this grain that is reduced to
   `T_module_out_sum` / `T_module_in_sum` / `T_network_inter_sum`. A companion
   long-format file would be needed to keep each ordered pair.
3. `Data_export_for_sharing/pca_metrics/pca_metrics_by_folder.csv` was **not**
   used: it covers only 2-module networks and its `pc1_var` / `n_genes`
   disagree with the `*_PCA_variance.csv` files, so it is an older lineage.
4. `_noPS` folders were not read (redundant — strict subset; use `has_PS`).
   `Single_EMT_MechSens` has no correlations output and is excluded.
