# Figure plan

Concrete, panel-by-panel plan for the interacting-teams paper. It takes the
five-figure skeleton from
[outline.md](Tkl_component_scatterplots/StateBased_newnorm/lmax10/outline.md)
and maps every panel onto the data file and script that can make it.

**Conventions fixed for the main text** (decided 2025-09-21):

| choice | value | why |
|---|---|---|
| path-length cutoff | `l_max = 10` | the endpoint of the published sweep; the `l_max` sweep goes to supplementary |
| phenotypic score | **influence-weighted** (`P_s_infl_*`) | r = 0.93-0.95 against network `T_s` with a null spread of ~0.045, against r ~ 0 and a null of ~0.33 for the state-based variant (`Metrics/loo_results.md`) |
| state aggregation | `top60`, with `top30` in supplementary | `dominant` is a sign-flip coin toss in 5 networks (near-tied mirror attractors, `Metrics/FORMULAE.md` §6) |
| **formalism** | **RACIPE leads Figure 1 (state heatmaps + PCA); Ising for the P_s panels downstream** | **revised 2026-09-24** — RACIPE is less ambiguous and carries fewer metrics, so the argument now opens with it. Supersedes the 2026-09-22 line below for Fig 1A/B/E and the supplementary module row. Previously: Every `P_s`, basin weight and state count comes from `<net>_finFlagFreq.csv`. The PCA metrics (`PC1_var`, `NumPC_90`, `PCA_entropy`, `PCA_entropy_norm`) come from `<Tag>_PCA_variance.csv`, which is RACIPE sampling (58k-265k samples/network) and has no Ising equivalent |
| corpus | the 56 two-to-five-module networks over AS, C, E, M, PS, RS | the only modules with a two-team split in `allNodes.nodes` |
| P_s sign | EMT oriented mesenchymal-positive; `TEAM1` in `Metrics/metrics_lib.R` | reproduces the published signs |

Any panel that departs from these is marked **[dev]** with the reason.

**Status legend**

| flag | meaning |
|---|---|
| **EXISTS** | a publication-usable figure file is already on disk |
| **REPLOT** | the data exist and are validated; only the plotting code is missing or needs restyling |
| **COMPUTE** | inputs exist but the quantity has never been computed |
| **MISSING** | the underlying simulation or network set does not exist yet |

---

## Asset inventory

What is already on disk, and what each thing is good for.

### Validated data tables

| file | grain | holds |
|---|---|---|
| `Metrics/metrics_network_level.csv` | 56 networks x 4 `l_max` | `T_network`, `T_inter_sum`, `T_network_intra_sum`, PCA metrics, topology metrics, all six network-level `P_s` columns |
| `Metrics/metrics_module_level.csv` | 56 x module x 4 `l_max` | `T_module_intra`, `T_module_out_sum`/`in_sum`, all six module-level `P_s` columns |
| `Metrics/metrics_pair_level.csv` | 56 x ordered module pair x 4 `l_max` | `T_from_to`, `T_to_from`, `T_inter_sum`, `T_intra_sum`/`mean` |
| `Metrics/module_properties.csv` | 6 modules x 4 `l_max` | `T_s_module`, `density_within`, `impurity_within`, intra/inter-team density, `density_across`, `impurity_across` |
| `Combined_data/combined_module_level_data.csv` | 1440 rows | the published-lineage module data, both `P_s` methods, with PCA joined |
| `Metrics/loo_correlations.csv` + `loo_null.csv` | removal x measure x `l_max` | leave-one-module-out Δr against a 2000-draw resampling null |
| `Metrics/clustering_inconsistency.csv` + `clustering_null*.csv` | module x k x `l_max` | three readings of clustering inconsistency, with nulls |

Every metric in the first four files is checked against the published
`correlations_*` outputs — the reproduction table is in `Metrics/README.md`.
The three documented departures from the formula image (abs-after-mean in
`T_kl`; node pairs not edges in `|E_kl|`; `N(N-1)` density) must be stated in
the methods.

### Simulation coverage

| formalism | role | coverage | location |
|---|---|---|---|
| **Ising** | **the formalism. Every phenotypic score, basin weight and state count.** | 28 two-module, 35 triple, 35 four, 21 five | `Two_module_nets/`, `Triple_EMT_MechSens/`, `Four_EMT_MechSens/`, `Five_EMT_MechSens/` (`*_finFlagFreq.csv`) |
| **RACIPE** | **the PCA spectrum, and nothing else** | all 56 networks | `*_EMT_Mec?Sens/<Tag>_PCA_variance.csv`; raw solutions in `6_module_racipe/` |
| logical rules | **out of scope** | rule files present for triples, fours, fives | `*_rules.txt` — never simulated |

The PCA metrics are the one place a second formalism enters, and they have no
Ising equivalent: `<Tag>_PCA_variance.csv` carries 58k-265k continuous samples
per network, against tens of discrete Ising attractors. Every panel using
`PC1_var`, `NumPC_90`, `PCA_entropy` or `PCA_entropy_norm` must say RACIPE on
its face; every panel using `P_s` must say Ising. `Figures/tier1_panels.R`
carries this as the `SRC_ISING` / `SRC_RACIPE` / `SRC_TOPO` caption constants.

Note also: `Boolean/module_network_dominant_results.csv` and
`6_module_racipe/module_ts_vs_pc1_analysis.csv` carry *identical* `basin_size`
values (Single_AS = 0.2935 in both), so **both are keyed to the same Ising
run** — the RACIPE file contributes `pc1_variance` only. Neither is an
independent formalism, and no cross-formalism P_s comparison exists.

Single-module Ising results appear in those two aggregate files
(`network_name = Single_AS`, ...) but the raw `Single_*_finFlagFreq.csv` are
not in the tree. They need to be recovered before Figure 1b can be built.

### Existing figure files

| location | n | content | reusable? |
|---|---|---|---|
| `Metrics/figures/` | 26 | LOO violins/impacts/properties, clustering readings A-C | yes, as Fig 4 source panels |
| `Cross_Module_Comparison/influence_heatmaps/` | 10 | module-grouped, team-segregated influence matrices, 3/4/5-module, `l_max` 10 | yes, directly — this is the Fig 2a style |
| `Cross_Module_Comparison/metric_clustering/` | 14 | clustermaps + dendrograms by module count | supplementary |
| `Cross_Module_Comparison/module_count_correlation/` | 11 | correlation vs number of modules, per metric | Fig 3 candidate |
| `Tkl_component_scatterplots/` | 108 | `T` components vs landscape outcomes | Fig 4 source, needs compositing |
| `entropy*/`, `scatter_all_params_*/` | ~430 | correlation heatmaps and the full scatter grid | supplementary / source |
| `panel1_plots/`, `panel2_plots/` | 33 | per-network `T_s` components vs `l_max` | diagnostic only — one exemplar at most |

---

## Figure 1 — Two-team modules and their phenotypic landscapes

*Message: each module on its own is a two-team system, and the strength of that
team structure sets how canalised its landscape is. This reiterates the known
result on this paper's modules, so it must be compact.*

### 1A — Module structure

| panel | content | status | source |
|---|---|---|---|
| **1a** | Influence matrices of one strong-team and one weak-team module, team-ordered, `l_max` 10 | **BUILT** `fig1a_influence_RS_PS.png` | `Single_EMT_MechSens/Single_{RS,PS}.topo` + `influence_full()`. RS shows clean +/- team blocks; PS is washed out |
| **1b** | Heatmap: 6 modules x properties (`T_s_module`, `density_within`, intra/inter-team density, `density_across`, `impurity_across`), column-scaled | **BUILT** `fig1b_module_properties.png` | `Metrics/module_properties.csv`, `l_max == 10` |

Two notes on 1a/1b:

- outline.md names **AS vs PS** as the strong/weak contrast. At `l_max = 10`
  the strongest module is actually **RS (0.917)**, then AS (0.766); PS is
  0.086. Use **RS vs PS** for the sharpest contrast, or keep AS vs PS and say
  why. Either way the caption should give both numbers.
- **Drop `impurity_within` from the 1b heatmap.** It is exactly 0 for five of
  six modules (PS alone at 0.143), so as a heatmap row it is a blank strip. It
  is worth one sentence of text — PS is the only internally impure module — not
  a row. Keep it in the supplementary version of the table.

### 1B — Landscape properties

| panel | content | status | blocker |
|---|---|---|---|
| **1c** | Steady-state heatmap (states x nodes, team-ordered) with an adjacent P_s strip | **BUILT** `fig1c_steady_states.png` | the **`Comb_<M>_<M>` self-pairs** in `Two_module_nets/` *are* the single-module runs — see below |
| **1d** | Per-state P_s by module, x ordered by `T_s_module`, **Ising only** | **BUILT** `fig1d_state_Ps_by_module.png` | same self-pair source. Drawn as basin-sized points, not a violin: each module has only 2-16 states |
| **1k** | **[dev]** isolated vs embedded commitment | **BUILT** `fig1k_isolated_vs_embedded.png` | not in outline.md; added because 1d showed every module commits alone, so the paper's question is what combination costs |
| **1e** | Frustration vs team strength | **BUILT** `fig1e_frustration.png` | **`frust0` was already a column of every `finFlagFreq.csv`** — written by the Boolean run. Nothing new computed. **r = -0.816**, the cleanest trend in the paper |
| **1f** | Hybrid-state frequency `F_hy` by module | **BUILT** `fig1f_hybrid_by_module.png` | computed from the per-state `P_s` + basin sizes. **Negative result — see below** |
| **1g** | Same violin, mean perturbation to switch `P_mean` on y | **BLOCKED — cluster** | needs a perturbation simulation that has never been run. `funcsKishore::retProbPert` exists as a helper but there is no pipeline |

### 1C — Trends against team strength

| panel | content | status | source |
|---|---|---|---|
| **1h** | `T_s_module` (x) vs P_s (y), facetted by `dominant` / `top60` | **BUILT** `fig1h_Ts_vs_Ps.png` | `Metrics/metrics_module_level.csv` |
| **1i** | team strength vs PC1 variance | **BUILT** `fig1i_Ts_vs_PC1.png` | `metrics_network_level.csv` — **[dev]** network level, not module (see below) |
| **1j** | `T_s_module` vs `F_hy` | **BUILT** `fig1j_hybrid_vs_Ts.png` | **Negative result — see below** |

**Recommended Figure 1 layout:** 1a | 1b on the top row, 1c across the middle,
1d | 1h | 1i on the bottom, with 1e promoted into the main figure — at
r = -0.816 it is the strongest single trend in the paper and it cost nothing.
1f/1j go to supplementary as a methods caveat, and 1g is cut unless the
perturbation run is worth the cluster time.

The earlier version of this plan called 1e-1g and 1j "three separate new
analyses". That was wrong for two of them: frustration was already computed
and sitting in the data, and hybrid states needed only arithmetic on scores
that already existed. Only 1g is genuinely new work.

---

## Figure 2 — Two-module combinations change the landscape

*Message: connecting two modules is not neutral. The same module pair can be
effectively uncoupled, coupled in a team-aligned way, or coupled in a
conflicting way, and the landscape follows.*

| panel | content | status | source |
|---|---|---|---|
| **2a** | Three influence matrices side by side: no interaction, team-aligned interaction, conflicting interaction | **REPLOT** | the plotting style already exists in `Cross_Module_Comparison/influence_heatmaps/` but only for 3-5 module networks; extend to 2-module |
| **2b** | 6x6 module-by-module heatmap, tile = `P_s_infl_dominant`, axes ordered by `T_s_module` | **BUILT** `fig2_pair_heatmaps.png` | `Metrics/metrics_network_level.csv`, `n_modules == 2`, `l_max == 10` |
| **2c** | Same, tile = `P_s_infl_top60` | **BUILT** | same |
| **2d** | Same, tile = `PC1_var` | **BUILT** | same |
| **2e** | Same, tile = `F_hy` | **BUILT as `PCA_entropy_norm`** | `F_hy` does not exist; `PCA_entropy_norm` carries the same "how concentrated is the landscape" meaning |

**Exemplar selection for 2a** — pick by data, not by eye. From
`metrics_pair_level.csv` at `l_max = 10`:

- *no interaction*: `Comb_AS_C` is the canonical case — exactly one
  cross-module edge (`Casp3 -| N_bcatenin`), `T_C->AS = 0`.
- *aligned* and *conflicting*: rank the 15 pairs by `impurity_across` and take
  the extremes. C (0.957) and M (0.842) are the misaligned connectors; AS
  (0.368), PS (0.357) and RS (0.360) the aligned ones. A concrete query over
  the pair table should pick the final two.

Note 2b-2e only use the 15 two-module networks, so each heatmap is a
15-cell upper triangle. Consider merging 2b-2d into one faceted panel.

---

## Figure 3 — Characterising inter-module interactions

*Message: here are the candidate metrics for "how strongly and how coherently
do two modules interact", and here is how they are distributed across the
corpus.*

| panel | content | status | source |
|---|---|---|---|
| **3a** | Violin: number of modules (x, 2-5) vs normalised cross-module edge density | **REPLOT** | `inter_density`/`density` in `metrics_network_level.csv`; the normalisation `n * density_across / sum(density_within)` still needs to be written |
| **3b** | Violin: n modules vs cross-module sign-coherence (the probability that edges from one team land with opposing signs) | **COMPUTE** | this is `impurity_across` generalised to the network level; `module_properties.csv` has it per module, not per network |
| **3c** | Violin: n modules vs `T_inter_sum` (influence-matrix inter-module strength) | **REPLOT** | `metrics_network_level.csv` |
| **3d** | Directed asymmetry `\|T_A->B - T_B->A\| / (T_A->B + T_B->A)` | **BUILT** `fig3d_asymmetry.png` | `metrics_pair_level.csv`. **New result:** asymmetry falls monotonically with module count, mean 0.79 (2 modules) -> 0.55 (5) |

`Cross_Module_Comparison/module_count_correlation/` already holds 11 plots of
correlation-vs-module-count per metric; those are the supplementary companion
to this figure.

**Worth testing while building 3b:** outline.md predicts the interaction-matrix
sign-coherence will be *less* informative than the influence-matrix `T_s`, and
suggests comparing the two for single modules as well. That comparison is a
single scatter (`impurity_across` vs `T_s_module`, n = 6) and belongs in this
figure as an inset.

---

## Figure 4 — Do the interaction metrics predict the landscape?

*Message: this is the paper's test. Score each interaction metric by how well
it predicts the emergent landscape.*

This is the figure with the most finished material — `Metrics/loo_results.md`
and `Metrics/clustering_results.md` are complete analyses with nulls.

| panel | content | status | source |
|---|---|---|---|
| **4a** | Baseline: `T_network` vs network `P_s_infl_top60`, all 56 networks, coloured by module count | **BUILT** `fig4a_Tnetwork_anchor.png` | `metrics_network_level.csv`; r = +0.943 |
| **4b** | `T_inter_sum` vs the same P_s, showing the *negative* relationship | **BUILT** `fig4b_Tinter_anchor.png` | same; r = -0.211 at l_max 10 |
| **4c** | r of 4b recomputed at each `l_max`: -0.54, -0.35, -0.29, -0.21 | **BUILT** `fig4c_lmax_trend.png` | **[dev]** new panel — at l_max 10 this correlation is at its weakest, so the sweep has to be shown |
| **4c** | Metric-scoring summary: each interaction metric vs each landscape outcome, tile = r, `l_max` 10 | **REPLOT** | `entropy_newnorm_stateps/correlations_pearson_equal_lmax10.csv` or recomputed from `metrics_network_level.csv` |
| **4d** | LOO: Δr by module against the resampling null band | **BUILT** `fig4d_loo_impact.png` | rebuilt at final size from `loo_correlations.csv` + `loo_null.csv` |
| **4e** | LOO violins: null distribution, six module removals, original r | **BUILT** `fig4e_loo_violins.png` | rebuilt at final size |
| **4f** | Clustering reading C: ARI of clustering vs contains-module, by `l_max` | **BUILT** `fig4f_clustering_ARI.png` | rebuilt at final size from `clustering_inconsistency.csv` |

**The honest finding this figure has to carry.** Two results are already
established and both are negative, so the figure must be designed to show them
rather than around them:

1. Against `T_network`, no module's removal moves the correlation beyond the
   resampling null (1 of 48 for the network-level P_s). The correlation is a
   property of the corpus, not of any module.
2. No module property predicts the LOO impact — every mean Spearman |ρ| is
   below 0.2, and `impurity_across` flips sign between predictors. With n = 6
   modules, the design cannot answer Q2.

The two things that *are* real: **M against the inter-module `T_s` sum**
(Δr up to 0.385, 8 of 48 removals clear the null) and **PS dominating the
clustering** (reading C, ARI 0.45 at `l_max` 10 against ~0 for all five
others, strengthening monotonically with `l_max`). Build 4d-4f around those
two and state the nulls in the caption.

`clustering_results.md` also flags a feature-set caveat that must be resolved
before this figure is final: the eight P_s columns are far from independent, so
the clustering is driven mostly by the influence-weighted pair. Rerun the
clustering on the two influence-weighted network-level columns only and confirm
PS still leads. **Do this before drawing 4f.**

---

## Figure 5 — Validation on artificial networks

*Message: build synthetic networks that match the biological ones on each
interaction metric, and check the metric still predicts the landscape.*

**Status: MISSING in full.** No artificial-network generator, network set, or
simulation exists in the repo. Every panel here is new work, and it depends on
Figure 4 having identified a metric worth validating.

Minimum viable version:

1. Generator: sample two-team modules with a target `T_s`, wire them with a
   target cross-module density and sign-coherence. The metric code in
   `Metrics/metrics_lib.R` can score candidates, so this can be
   rejection-sampled.
2. ~50 networks per condition across a 2-D grid of (cross-module density,
   sign-coherence), Ising only.
3. Panels: 5a the design grid; 5b landscape outcome across the grid; 5c the
   biological 56 overlaid on that grid; 5d predicted-vs-observed P_s.

outline.md already predicts the proposed metrics will *not* predict well.
Figure 4's results support that prediction. **Consider whether Figure 5 is
worth building before Figure 4's metric scoring has identified something to
validate** — if it has not, Figure 5 becomes a large simulation effort to
confirm a negative, and the paper may be better served by the Figure 3/4 pair
plus an expanded corpus.

---

## Supplementary

| S | content | status |
|---|---|---|
| S1 | Module property table including `impurity_within`, all `l_max` | REPLOT |
| S2 | `l_max` sweep: every main-text correlation at `l_max` 1, 5, 7, 10 | EXISTS (`entropy*/heatmaps/`) |
| S3 | State-based P_s versions of Figs 2 and 4, with the null spreads | REPLOT — the variant disagreement is itself a result |
| S4 | Full scatter grid, all metric-outcome pairs | EXISTS (`scatter_all_params_*/`, 240 files) |
| S5 | `T_s` component decomposition vs `l_max`, one exemplar network | EXISTS (`panel1_plots/`) |
| S6 | Clustering dendrograms and clustermaps by module count | EXISTS (`Cross_Module_Comparison/metric_clustering/`) |
| S7 | LOO and clustering panels not used in Fig 4 (20 of the 26) | EXISTS (`Metrics/figures/`) |
| S8 | Methods: the three departures from the formula image, and the reproduction table | EXISTS (`Metrics/FORMULAE.md`, `Metrics/README.md`) |
| S9 | Methods: why PCA is RACIPE-derived while everything else is Ising, and why the two are not interchangeable | TEXT — no figure needed |

---

## Work queue

Ordered by ratio of paper value to effort.

**Tier 1 — DONE.** `Rscript Figures/tier1_panels.R` writes all of it to
`Figures/panels/`: eleven panels (1b, 1h, 1i, 2b-2e, 3d, 4a-4f) plus
`FIG4_composite.png` and `FIG1_partial_composite.png`. Departures and new
findings are recorded in the section below.

**Tier 2 — needs new computation, inputs exist**

7. Fig 3b: network-level cross-module sign-coherence.
8. Fig 3a: the normalised cross-module density.
9. Fig 2a: extend the influence-heatmap code to two-module networks, and pick
   the three exemplars from `metrics_pair_level.csv` by query. The
   single-module version in `Figures/fig1_remaining.R` is the starting point.
10. Fig 4f prerequisite: rerun clustering on the influence-weighted columns
    only, confirm PS still leads.

**Tier 3 — one item left**

12. Perturbation analysis for `P_mean` → unblocks 1g, if wanted at all. This
    is the only panel in the paper still needing a simulation that does not
    exist.

~~13. RACIPE P_s pipeline over `6_module_racipe/*_solution*.dat`.~~
~~14. Run the logical formalism from the `*_rules.txt` files.~~
**Both dropped 2026-09-22** by the Ising-only decision. RACIPE stays in the
paper for the PCA spectrum, which it alone provides.

**Tier 4 — new analyses, scope decision required**

15. Frustration (1e), hybrid states (1f, 1j, 2e), perturbation (1g).
16. The whole of Figure 5.

---

---

## Build record — Tier 1

`Figures/tier1_panels.R` -> `Figures/panels/`. Re-runnable from the repo root;
it reads only the validated `Metrics/*.csv` tables.

### Departures from the plan above

1. **1i is plotted at the network level, not the module level.** `PC1_var` is a
   property of the network, so joining it onto the 180 module rows would repeat
   each y value three to five times over and inflate n. 56 clean points,
   r = +0.614.
2. **2e uses `PCA_entropy_norm` in place of `F_hy`**, which does not exist.
3. **4c is a new panel.** At `l_max` 10 — the convention everywhere else — the
   `T_inter_sum` correlation is at its *weakest* (-0.21 against -0.54 at
   `l_max` 1). Showing only the `l_max` 10 scatter would understate the result,
   so the sweep gets its own panel. It was first drawn as an inset inside 4b
   and collapsed illegibly; insets on a dense scatter do not survive.
4. **4d uses the pooled 95th percentile of |null deviation| (0.210)**, the same
   criterion as `loo_results.md`, rather than a per-`l_max` central-90% band.
   The per-`l_max` band is tighter at low `l_max` and would have flagged 12 of
   48 removals against the doc's 7.

### Corrections the data forced during the build

- **4d's title was wrong on the first pass.** "Only M's removal moves this
  correlation beyond chance" does not hold: 7 of 48 removals clear the null —
  M in four cases (max Δr = +0.328) and **AS in three** (max -0.211). The
  panel now names both.
- **3d's caption was wrong on the first pass.** 4 of 225 unordered module pairs
  have no coupling in either direction and are dropped, not 8 — the pair table
  stores both orderings, so the raw NA count double-counts. The count is now
  computed in the script rather than written by hand.
- **1h was misleading as first drawn.** Module team strength is clustered into
  six module-specific bands on x, so a single fit over 180 points reads as far
  more evidence than there is — and the slope is largely carried by PS sitting
  alone at 0.09. The panel now colours and direct-labels by module, marks the
  six module means, and reports the mean *within-module* r (+0.39 dominant,
  +0.45 top60) beside the overall r (+0.652, +0.714). The within-module
  correlation is the honest number, and it holds up.

### Formalism labelling (2026-09-22)

Every panel now states its source on its face: `SRC_ISING` on the P_s panels
(1h, 2b, 2c, 4a-4f), `SRC_RACIPE` on the PCA panels (1i, 2d, 2e), `SRC_TOPO` on
the purely structural ones (1b, 3d). Figure 2 mixes both, so its panel titles
carry the formalism inline — "P_s, top 60% (Ising)", "PC1 variance (RACIPE)".
No numbers changed: the panels already drew P_s from Ising and PCA from RACIPE.

### New result

**Inter-module coupling is strongly one-directional, and less so as networks
grow** (3d). Mean directed asymmetry falls monotonically with module count:
0.79 (2 modules), 0.73 (3), 0.64 (4), 0.55 (5), where 1 is strictly one-way and
0 perfectly reciprocal. This is the "any other metrics?" slot in outline.md and
it is the one interaction metric here with a clean monotone trend, so it is
worth carrying into the Fig 4 metric scoring.


---

## Build record — Figure 1 remainder

`Figures/fig1_remaining.R` -> `Figures/panels/`. Four more panels, none of
which needed a new simulation.

### Two panels I had wrongly marked MISSING

1. **Frustration (1e) was already in the data.** Every `finFlagFreq.csv`
   carries a `frust0` column written by the Boolean run — populated in 117 of
   119 networks, range 0 to 0.405. No new computation was needed at all.
   Basin-weighted mean frustration against network team strength gives
   **r = -0.816 over the 56 networks**, the cleanest trend anywhere in this
   project and stronger than anything in Figure 3 or 4. It should be promoted
   into main-text Figure 1.
2. **Hybrid states (1f, 1j) needed only arithmetic**, not simulation: the
   per-state `P_s` and basin sizes were already there via
   `metrics_lib.R::state_ps_table()`.

### Negative result — `F_hy` as defined does not work

Defining a hybrid state as `|P_s^state| < threshold` fails two ways:

| threshold | r(F_hy, team strength) | r(F_hy, module size) |
|---|---|---|
| \|P_s\| < 0.5 | -0.22 | -0.09 |
| \|P_s\| < 1.0 | +0.29 | **-0.88** |
| \|P_s\| < 1.5 | -0.15 | -0.36 |

- **It does not track team strength.** RS, the *strongest* module (T_s 0.92),
  is the most hybrid at every threshold.
- **It tracks module size instead**, at r = -0.88 at the main threshold. A
  fixed cut on `P_s` is not size-invariant: a 6-node module (M) quantises
  `P_s` far more coarsely than a 17-node one (AS), so it falls inside the
  hybrid band far more often.
- **The module ranking is not stable across thresholds** — PS is 2nd at 0.5,
  3rd at 1.0, 2nd at 1.5; M is 5th at 0.5 and 2nd at 1.0.

outline.md anticipated exactly this ("hybrid states defined strictly on the
phenotypic score to start with. The definition can be revisited based on the
results"). It needs revisiting: a size-invariant definition — per-node
committed fraction, or a cut calibrated to each module's attainable `P_s`
lattice — before `F_hy` can carry any argument. Per-instance values for all
three thresholds are in `Figures/panels/fig1f_hybrid_sensitivity.csv`.

### 1a

RS against PS, rather than outline.md's AS against PS: RS (0.92) is the
strongest module and the block structure is unmistakable, while PS (0.09) is
visually washed out. AS (0.77) is second-strongest and the caption says so.

---

## The single-module runs — found, and where

**They are the `Comb_<M>_<M>` self-pairs in `Two_module_nets/`.** I first
reported them missing; that was wrong. `compute_metrics.R` deliberately drops
self-pairs ("drop self-pairs (single modules)"), and searching for
`Single_*finFlagFreq*` by name found nothing, so they were invisible to both
the pipeline and to me.

Verified 2026-09-22, all six modules: topology identical to `Single_<M>.topo`
edge for edge, and `Comb_<M>_<M>_nodes.txt` identical to `Single_<M>_nodes.txt`.

| module | nodes | edges | states |
|---|---|---|---|
| AS | 17 | 52 | 8 |
| C | 14 | 22 | 2 |
| E | 12 | 40 | 2 |
| M | 6 | 9 | 6 |
| PS | 15 | 68 | 16 |
| RS | 7 | 19 | 6 |

### Module definitions cross-checked

Against `allNodes.nodes` and the full topology
`../../mechanosensing/EMT_Mechanosensing_TGFbeta.topo` (633 edges, 150 nodes;
`allNodes.nodes` has exactly 150 rows). For all six modules the node set from
`allNodes.nodes` `Program`, from `Single_<M>_nodes.txt`, and from the
self-pair agree exactly, and every module node appears in the full topology.
Team splits: AS 13/4, C 8/6, E 9/3, M 5/1, PS 12/3, RS 4/3 — matching
`Metrics/README.md`.

### Still blocked

Only **1g (`P_mean`)**, which needs a perturbation simulation that has never
been run. `funcsKishore::retProbPert` exists as a helper but there is no
pipeline. Nothing else in Figure 1 is blocked.

Cluster state for the record (Explorer, `/scratch/a.hari/InteractingTeams/`):
the same 119 `finFlagFreq.csv` files as the local tree, plus `PCA_variance/`,
`run_boolean.sh` and `run_racipe.sh`. No extra simulations live there.

**A near-miss to avoid.** `../Singles_MiDAS/` on this mac also holds
`Single_*_finFlagFreq.csv` for all six modules — but it is a different
project's module set. Only C and PS have identical topologies; AS, E, M and RS
differ. Its `Single_M` is not a bigger Migration module but a **different
module entirely** — 20 nodes of mitochondrial dynamics (`NADp_m`,
`TCA_cycle`, `ETC_Normal`, `Drp1`, `MFN1_2`, `PINK1`, `mROS`), against
Migration's 6 nodes (`PAK1`, `Merlin`, `Rac1`, `IQGAP1_LeadingE`,
`Horizontal_Pol`, `Fast_Migration`). Same letter, unrelated network.
**Do not use it for this paper.**


---

## Build record — single-module panels

`Figures/fig1_singles.R` -> 1c, 1d, 1k. Source: the `Comb_<M>_<M>` self-pairs.

### The headline this produced

**In isolation, every module commits — team strength makes almost no
difference.** Basin-weighted mean |P_s| spans only 1.72 to 2.00 across the six
modules and correlates -0.10 with team strength. C and E are perfectly
bistable (2 states, |P_s| = 2 in 100% of the basin); even PS, the weakest
module at T_s 0.086, sits at 1.83.

**What team strength predicts is how much a module loses when it is
combined** (panel 1k):

| module | T_s | isolated | embedded (mean of 30) | loss |
|---|---|---|---|---|
| PS | 0.086 | 1.83 | 1.67 | **+0.16** |
| C | 0.497 | 2.00 | 1.93 | +0.07 |
| M | 0.576 | 1.74 | 1.84 | **-0.10** |
| E | 0.715 | 2.00 | 1.99 | +0.01 |
| AS | 0.766 | 1.91 | 1.81 | +0.10 |
| RS | 0.916 | 1.72 | 1.71 | +0.01 |

r(T_s, isolated) = **-0.10**, r(T_s, loss) = **-0.51**. Weak-team modules lose
the most; M is the exception and actually gains. With n = 6 this is
suggestive, not significant — but it is the paper's thesis stated at the
module level, and it motivates Figures 2-4 directly.

### A better hybrid metric than the one that failed

1c gives a size-invariant alternative to the `|P_s| < threshold` definition
that broke in 1f/1j: **the basin fraction in fully-committed states**
(|P_s| = 2 exactly).

| module | PS | AS | M | RS | C | E |
|---|---|---|---|---|---|---|
| committed basin fraction | **0.139** | 0.584 | 0.743 | 0.783 | 1.000 | 1.000 |

PS is an order of magnitude below every other module, and its *dominant* pair
is not even fully committed (|P_s| = 1.83 against 2.00 for all five others).
This is the cleanest separation PS shows anywhere, and it is consistent with
PS dominating the clustering in Figure 4f. Worth adopting in place of `F_hy`.

### Note on 1d's form

A violin is wrong here: the modules have 2, 2, 6, 6, 8 and 16 states. The
panel plots the states themselves, sized by basin, with the basin-weighted
mean |P_s| marked.


---

## Correlation heatmaps — inputs x outputs by module count

`Figures/correlation_heatmaps.R` -> `Figures/correlation_heatmaps/`. Two
figures plus their underlying `r`/`p` tables, at `l_max = 10`:

| file | outputs |
|---|---|
| `corr_heatmap_infl_lmax10.png` | influence-weighted P_s |
| `corr_heatmap_state_lmax10.png` | state-based P_s |
| `corr_heatmap_<variant>_lmax10.csv` | every `r`, `p` and `n` behind the tiles |

Eight structural inputs (network `T_s`, sum of module `T_s`, inter-module
`T_s` sum, impurity, density, intra/inter-team density, network size) against
eleven landscape outputs (network and module P_s at dominant/top30/top60, PC1
variance, PCs to 90%, PCA entropy, frustration, number of states), faceted by
module count.

**n = 1 is now a real column.** It comes from
`Metrics/metrics_*_singles.csv`, written by
`Rscript Metrics/compute_metrics.R 10 --singles` — a new opt-in flag that runs
the same validated driver over the `Comb_<M>_<M>` self-pairs and writes to
suffixed files, leaving the 56-network tables untouched. Group sizes are
6 / 15 / 20 / 15 / 6.

### What the two figures show

**The two P_s variants are not comparable as readouts.** Over the P_s output
rows, the influence-weighted figure has **102 of 234 cells significant at
p < 0.05; the state-based figure has 1**. This is the same split
`loo_results.md` found, now resolved by module count, and it is the clearest
argument for the main-text convention.

Note that five of the eleven output rows — PC1 variance, PCs to 90%, PCA
entropy, frustration, number of states — do not involve P_s and are therefore
**identical in both figures**. Only the six P_s rows differ.

### Two results that hold across every module count

1. **Frustration against the summed module team strength** strengthens
   monotonically with network size and is the most consistent relationship in
   the whole dataset:

   | n modules | 1 | 2 | 3 | 4 | 5 |
   |---|---|---|---|---|---|
   | r | -0.64 | **-0.86** | **-0.93** | **-0.96** | **-0.97** |

   It beats network `T_s` as a predictor of frustration at every n >= 2, and
   it is stronger than anything in Figures 3 or 4.

2. **Influence-weighted network P_s against network `T_s`** is flat and high
   at 0.95-0.998 for n = 1 to 4, then drops to 0.75 at n = 5 (p = 0.09, n = 6
   — so the drop is not established, but it is the one place the anchor
   relationship weakens).

### Caveats written onto the figures

- **The n = 1 network `T_s` / network P_s cell (r = 0.998) is close to an
  identity, not a finding.** For a fully committed state P_s = 2 x `T_s`, and
  in isolation almost every module is fully committed (see the single-module
  build record). The observed ratio is 2.006 for RS at n = 1 but spans
  1.16-2.04 once two modules are combined — so the relation is tautological
  in isolation and only becomes informative under combination.
- **`T_inter_sum` is identically 0 at n = 1**, so that column is blank (grey)
  rather than zero-correlation.
- **Sample sizes are small.** At n = 6 a correlation needs |r| > 0.81 for
  p < 0.05; at n = 20, |r| > 0.44. Unstarred cells are noise. The n = 1 and
  n = 5 panels have six networks each and should not be read as trends.


### Do they match the published correlation tables?

Checked with `Figures/correlation_heatmaps/verify_against_published.R`, which
compares every reproducible cell against
`*/correlations_{double,triple,four,five}{,_stateps}/correlations_pearson_equal_lmax10*.csv`.
504 cells. **Yes for influence-weighted, with two documented exceptions; no
for state-based, by design.**

**Influence-weighted — exact where it should be:**

| output row | cells | median \|diff\| | worst | within 0.01 |
|---|---|---|---|---|
| PC1 variance | 28 | **0** | 0 | 28/28 |
| PCs to 90% | 28 | **0** | 0 | 28/28 |
| PCA entropy (norm) | 28 | **0** | 0 | 28/28 |
| module P_s, dominant | 28 | **0** | 0 | 28/28 |
| network P_s, dominant | 28 | **0** | 0 | 28/28 |
| module P_s, top30 | 28 | 0.014 | 0.093 | 8/28 |
| network P_s, top30 | 28 | 0.041 | 0.094 | 4/28 |
| module P_s, top60 | 28 | 0.014 | 0.053 | 12/28 |
| network P_s, top60 | 28 | 0.017 | 0.051 | 11/28 |

Everything that does not depend on basin weights reproduces to zero. The
top30/top60 rows are off by at most 0.094 — this is the drift already
recorded in `Metrics/FORMULAE.md` §8: the `Avg0` values in the current
`finFlagFreq.csv` are not the ones the published columns were computed from
(ordinary Monte-Carlo scatter between Boolean runs). It is a property of the
inputs, not of this code.

**State-based — different by construction, and the difference is understood.**
The published `_stateps` columns carry an extra `1/N_network` factor
(`Metrics/README.md`, known gap 2). Because `N` varies across networks that
is a per-point rescaling, so unlike a global constant it **does** change a
correlation. Dividing my scores by `n_nodes` recovers the published values for
the `dominant` rows:

| row | median \|diff\| as-is | median \|diff\| after /N | within 0.01 after /N |
|---|---|---|---|
| module P_s, dominant | 0.032 | **0.000** | 21/28 |
| network P_s, dominant | 0.030 | 0.004 | 19/28 |
| module P_s, top30 | 0.225 | 0.219 | 0/28 |
| network P_s, top30 | 0.255 | 0.272 | 0/28 |
| module P_s, top60 | 0.217 | 0.208 | 1/28 |
| network P_s, top60 | 0.208 | 0.210 | 2/28 |

So the `dominant` rows agree once the factor is accounted for, and the
top30/top60 rows do not reproduce under any constant factor — exactly what
`Metrics/README.md` gaps 2 and 3 predict, since the old aggregates were
computed differently *and* on drifted basin weights.

**The state-based heatmap therefore deliberately does not reproduce the
`_stateps` folders.** It follows the stated definition, with the `1/N` dropped
per the instruction recorded in `FORMULAE.md` §6. If the published lineage is
wanted for comparison, divide by `n_nodes`. This is a live choice, not a
settled one — see open decision 6.


---

## The two P_s definitions, and a caveat on Figure 4a

Asked 2026-09-23; both answers are now in `Metrics/FORMULAE.md` §6 and §7.

### State-based (§6) — range [-2, 2], and why

`decode_state()` codes each node `+/-1`, so a team mean is `2*E(T) - 1` where
`E(T)` is the activity fraction the formula sheet uses. The difference is
therefore **exactly twice** the sheet's score:
`mean(S over T1) - mean(S over T2) = 2*(E(T1) - E(T2))`. The sheet's range is
[-1, 1]; the implementation's is [-2, 2]. This is a **fourth departure** from
the formula image, now documented (the other three are in §3 and §5).

Correlations are invariant to the factor of 2, so nothing in
`Metrics/*.csv` or `Figures/correlation_heatmaps/` changes. Absolute cuts do
change: "fully committed" is 2 here and 1 on the sheet, and the hybrid
thresholds 0.5/1.0/1.5 are 0.25/0.5/0.75 in sheet units. The affected
captions now say so.

### Influence-weighted (§7) — range [0, 2], observed 0.01-1.62

`F_T` is the mean over the team of each node's outgoing influence, weighted by
the target's state and averaged over all `N` network nodes;
`P_s^infl = |F_T1 - F_T2|`. It is a *magnitude* — the absolute value discards
the direction the state-based score keeps.

### The caveat this surfaced — Figure 4a is largely structural

`P_s^infl` multiplies the state by `I`, and on a fully committed state it
reduces to a weighted sum of the same four block means that define
`T_network`. Measured:

| quantity | r vs `T_network` (56 networks, l_max 10) |
|---|---|
| `P_s^infl` on a **synthetic** committed state — no simulation | **0.996** |
| `P_s^infl` on the actual dominant state (this is Fig 4a) | 0.934 |
| `P_s^state` on the actual dominant state | 0.044 |

In **30 of 56 networks the dominant state is the committed state**, so there
the simulated and synthetic numbers are identical. On the other 26 the actual
correlation drops to **0.564** while the synthetic holds at 0.975.

**Figure 4a's r = 0.94 is therefore not an emergent simulation result — it is
mostly an algebraic consequence of both axes being built from the same
influence matrix.** This also explains the n = 1 near-identity already flagged
on the correlation heatmaps, and it reframes the influence-vs-state gap: the
state-based score's r ~ 0 is not a defect of that measure, it is what the
simulation actually says once the structural term is removed.

### How much survives — the partial correlation

Controlling for the synthetic-committed baseline:

| | r vs `T_network` | partial r, controlling for synthetic |
|---|---|---|
| `P_s^infl` (dominant) | +0.934 | **+0.148** |
| `P_s^state` (dominant) | +0.044 | -0.121 |

**Essentially nothing survives.** The r = 0.93 is the structural term.

### The one genuinely dynamical result inside it

Whether the dominant state *is* the committed state is not algebra — it is an
outcome of the Boolean dynamics, and team strength predicts it:

| | committed dominant state | not committed |
|---|---|---|
| networks | 30 | 26 |
| mean `T_network` | **0.474** | **0.284** |

point-biserial r = **+0.570**, p < 1e-4. *This* is the defensible version of
"team strength shapes the landscape": stronger teams make the committed state
win, rather than making a continuous score larger.

### Why widening the state window changes r — same cause

`dominant` -> `top30` -> `top60` does not move r in one direction; it depends
on whether the dominant state was already committed:

| n modules | committed dominant | r(dominant, synthetic) | r(top60, synthetic) |
|---|---|---|---|
| 2 | 9/15 | 0.999 | 0.949 |
| 3 | 12/20 | 0.942 | 0.955 |
| 4 | 7/15 | 0.861 | **0.961** |
| 5 | 2/6 | **0.225** | **0.779** |

As module count rises the dominant state is less and less often the committed
one, and `r(dominant, synthetic)` collapses from 0.999 to 0.225. But
`r(top60, synthetic)` barely moves. The committed mirror pair still holds most
of the basin even when neither state is strictly top, so averaging over the
top 60% pulls the estimate back onto the committed configuration.

So widening the window **is** the same circularity, recovered rather than
diluted: it does not add independent information, it makes the estimate of the
*same structural quantity* less noisy. That is why `T_s` vs network
`P_s^infl` rises from 0.206 (dominant) to 0.747 (top60) at n = 5, and falls
slightly at n = 2 where the dominant state was already committed.

**One exception worth keeping separate:** the *state-based* module-level
correlation rises with the window at every module count (n = 1 to 5, all up).
That score contains no influence matrix, so it cannot be the structural
effect — there, widening is genuine noise reduction on a near-saturated
measure (most attractors sit at |P_s| = 2, so the dominant state alone has
little variance to correlate with).

**Actions this implies, not yet taken:**

1. Re-caption Fig 4a against the synthetic-committed baseline rather than
   against zero, or demote it from "the anchor result".
2. Re-run the Figure 4 metric scoring with the synthetic-committed `P_s` as a
   null, so `Delta r` is measured against the structural expectation.
3. Revisit whether the influence-weighted score should be the main-text
   convention at all. `loo_results.md` chose it for its tight null and high
   r — both of which are now partly explained by circularity.

This is a paper-level decision, so nothing has been changed on the basis of
it yet; see open decision 7.


---

## Interaction-weighted P_s — a third variant, and what it settles

Added 2026-09-23 on request. "Interaction matrix" is this codebase's name for
the **signed adjacency `A`** (`funcsKishore::TopoToIntMat`), so the
interaction-weighted score is §7 with `A` in place of `I`:

```
F_T = (1/|T|) * sum_{i in T} ( (1/N) * sum_j A[i,j] * S_j ),   P_s^int = |F_T1 - F_T2|
```

**It needed no new estimator.** `FORMULAE.md` §1 states `I = A` exactly at
`l_max = 1`; verified here to 0 difference on the matrices and to 6 decimals
on the score (`Comb_AS_C`/AS: 0.248759 both ways). So the interaction-weighted
column *is* the `l_max = 1` influence column, and
`Figures/correlation_heatmaps.R` reads it from the `l_max = 1` rows rather
than recomputing. It is **l_max-free** by construction — `A` has no path
cutoff — while the structural inputs stay at `l_max = 10`.

Third figure: `corr_heatmap_int_lmax10.png` (138 of 429 cells significant,
against 163 influence and 62 state-based).

### It does not fix the circularity

| P_s variant | r vs `T_network` (l_max 10) | its own synthetic-committed baseline | partial r |
|---|---|---|---|
| influence-weighted (`I`, l_max 10) | +0.934 | +0.996 | **+0.148** |
| **interaction-weighted (`A`)** | **+0.503** | **+0.626** | **-0.170** |
| state-based | +0.044 | — | -0.121 |

The raw correlation drops from 0.93 to 0.50, which looks like progress, but
the partial correlation against its *own* structural baseline is **-0.170** —
nothing survives there either. The baseline (0.626) is again *higher* than the
observed value (0.503).

The drop is a scale mismatch, not a fix: `A` and `I(l_max 10)` are different
matrices, so pairing an `A`-based score with an `I`-based `T_s` simply
weakens the shared term. Pair it with its own-scale team strength and the
circularity returns in full: **r(`T_s` at l_max 1, interaction P_s) = +0.916**.

**Conclusion — revised 2026-09-23.** "Circularity" was the wrong frame. The
matrix-weighted scores are **team-structure alignment measures**, and exactly
so: `P_s = |w^T M S| / N` with `w` the size-normalised team vector, verified
to 2.8e-17. Frustration is the same bilinear form with the state in the left
slot, `S^T A S` (also exact). Both project `M*S` — one onto the state, one
onto the team partition — and they coincide when the state *is* the partition,
which is why 30 of 56 networks sit on the synthetic-committed value.

So the adjacency/influence P_s is not a corrupted landscape readout; being a
joint function of topology and state is the point. It is finer than
frustration: team-resolved, signed, and per-team normalised.

What this changes: correlating it with `T_network` is a fair question — *do
attractors align with the team partition, more so when teams are stronger?*
— and the answer is yes. It does **not** license the looser claim that team
strength predicts a *landscape outcome*, since both sides derive from the same
topology (partial r = +0.148 against the synthetic baseline). Report it as an
alignment result. The state-based score remains the one to use when a readout
independent of topology is wanted, and its r ~ 0.04 is then a measurement
rather than a defect. See `Metrics/FORMULAE.md` §7.

All three heatmaps now carry this caution on their face, and the two
matrix-based ones name their partial correlations.


---

## Inter-module team strength — the definition was not team-aware

Raised 2026-09-24. `T_inter` (§3.4) averaged `I` over the **whole A x B
block**, ignoring the team split inside each module. Two separate defects
follow, and both bias the result.

### Defect 1 — cancellation

A coupling that is perfectly aligned with team structure (A's team 1 driving
B's team 1 and opposing B's team 2) has sub-blocks of opposite sign, which
cancel in a whole-block mean. **The most coherent coupling possible reports
near zero.**

Over the 450 ordered module pairs at `l_max` 10, the team-resolved value is
larger in **361 (80%)**, median **4.3x**, 90th percentile **29x**. The extreme
is M->C: whole-block **0.0714** against team-resolved **0.740**, with
`T_aligned = -0.740` — i.e. all four sub-blocks carry the aligned pattern at
full strength, a perfectly coherent (anti-aligned) coupling read as ~0. Old
and new correlate only **+0.61**, so this reorders pairs rather than rescaling
them.

M being the worst case is worth noting: M is also the module whose removal
dominated the leave-one-out analysis against `T_inter_sum`
(`loo_results.md` Q1). That effect may simply be M's coupling being
mismeasured.

### Defect 2 — the roll-up is a sum over n(n-1) pairs

`T_inter_sum` adds over every ordered pair, so it grows mechanically with
module count: **r(T_inter_sum, n_modules) = +0.81**. That confound inverts the
pooled correlation — a textbook Simpson's paradox.

### Both fixed: the headline reverses

r against network `P_s` (influence-weighted, top60), `l_max` 10:

| metric | pooled (56) | n=2 | n=3 | n=4 | n=5 | r with n_modules |
|---|---|---|---|---|---|---|
| `T_inter_sum` (current) | **-0.211** | +0.570 | +0.469 | +0.569 | +0.162 | +0.81 |
| current, mean-normalised | +0.403 | +0.570 | +0.469 | +0.569 | +0.162 | +0.10 |
| team-blocks, sum | -0.129 | +0.869 | +0.920 | +0.854 | +0.470 | +0.82 |
| **team-blocks, mean** | **+0.680** | **+0.869** | **+0.920** | **+0.854** | **+0.470** | **+0.11** |

**Within every module count the relationship was already strongly positive**;
the published negative pooled value is an artefact of the sum-over-pairs
confound, and the team-blindness was costing roughly 0.35-0.45 of correlation
on top. Corrected, inter-module coupling is a *strong positive* predictor of
the phenotypic score (0.85-0.92 within module count), not the weak negative
one reported in `loo_results.md`.

### What was implemented

`Metrics/metrics_lib.R` gains two functions, documented there in full. With
`m_kl = mean(I[A_k, B_l])` over the four team sub-blocks:

- `T_inter_blocks  = mean(|m_11|, |m_12|, |m_21|, |m_22|)` — magnitude of
  team-resolved coupling; the direct analogue of `T_network()`, which averages
  four `|blocks|` the same way.
- `T_inter_aligned = (m_11 - m_12 - m_21 + m_22)/4` — *signed* coherence with
  the team structure, and exactly `w_A^T I w_B / 4` for the size-normalised
  team vectors. This is the inter-module member of the same bilinear family as
  the phenotypic scores, and it is what outline.md's Figure 3b asked for
  ("the probability that edges from one team to another have opposite signs").

`Metrics/compute_inter_teamaware.R` writes
`metrics_inter_teamaware.csv` (1800 pair rows) and
`metrics_inter_teamaware_network.csv` (224 network rows, with sum, mean and a
`T_inter_coherence = sum|aligned| / sum blocks` ratio).

**The validated tables are untouched.** `T_inter` still reproduces the
published values 1800/1800 and stays as `T_inter_old` in the new file, so the
published lineage remains checkable.

`T_inter_coherence` is `NA` for `Comb_M_PS` and `Comb_M_RS` — verified to have
**zero** cross-module edges, so it is a genuine 0/0, not a bug.

### Not yet propagated

The correlation heatmaps, Figure 3 and Figure 4 still use the old
`T_inter_sum`. Rebuilding them on `T_inter_blocks_mean` is the obvious next
step, and it will change Figure 4's sign. The leave-one-out analysis should
be rerun on it too, since M's apparent dominance may be an artefact.


---

## Correlation matrices, reclassified (2026-09-24)

`Figures/correlation_heatmaps.R` now writes **six** matrices to
`Figures/correlation_heatmaps/`, classified by output family so each has one
grain and one readout. The previous three-file layout is deleted on rerun.

| # | file | outputs | grain |
|---|---|---|---|
| 1 | `corr_pc_lmax10` | PC1 var, PCs to 90%, PCA entropy (raw + norm), frustration, state count | network |
| 2 | `corr_module_infl_lmax10` | module P_s, influence — mean **and variance** over modules | network |
| 3 | `corr_module_state_lmax10` | module P_s, state — mean **and variance** | network |
| 4 | `corr_network_infl_lmax10` | network P_s, influence | network |
| 5 | `corr_network_state_lmax10` | network P_s, state | network |
| 6 | `corr_modulegrain_lmax10` | each module's own P_s (both variants) | **module** |

Matrix 1 holds the outputs that do not depend on the P_s definition, so it is
shared by every variant rather than repeated in each — that was the redundancy
in the old layout.

Inputs now include the **team-aware** `T_inter_blocks_mean` and the signed
`T_inter_aligned_mean`, with the old `T_inter_sum` retained as a labelled row
for continuity.

### Module-level variance — a clean negative

The variance of P_s across a network's modules was added to matrices 2 and 3
to ask whether structure predicts how *unevenly* the modules score, which the
mean cannot see. **It does not: essentially no cell reaches p < 0.05 in either
variant, at any module count.** Whatever sets the spread between modules
within a network is not in this set of structural predictors. (Variance is
`NA` at n = 1 and rests on two values at n = 2.)

### Matrix 6 — the module-grain result

Both sides module-level, so no roll-up and no repeated predictor, and the
samples are 2-10x larger: **6 / 30 / 60 / 60 / 30**.

**Outgoing inter-module coupling predicts a module's own phenotypic score
better than its own intra-team strength does, and incoming coupling does not
predict it at all.** Against module P_s (influence, top60):

| n | rows | r intra | r **out** | r in | partial out\|intra | partial intra\|out |
|---|---|---|---|---|---|---|
| 2 | 30 | +0.744 | **+0.770** | +0.323 | +0.822 | +0.803 |
| 3 | 60 | +0.715 | **+0.895** | +0.266 | +0.876 | +0.658 |
| 4 | 60 | +0.795 | **+0.912** | +0.247 | +0.854 | +0.640 |
| 5 | 30 | +0.866 | **+0.926** | +0.184 | +0.813 | +0.635 |

Both survive controlling for each other, and in a joint model both betas are
independently significant at every module count (p < 1e-4 throughout,
R² = 0.86-0.92). The in/out asymmetry is the substantive finding: **a module's
phenotype tracks how strongly it drives the rest of the network, not how
strongly it is driven.**

Three reasons this was invisible before: the network-grain roll-up averages it
away, the old `T_inter` cancelled team-aligned coupling toward zero, and
pooling across module counts hides it (pooled r_out = +0.468 with
partial +0.199, because collinearity between `out` and `intra` rises from 0.34
at n = 2 to 0.77 at n = 5).

As always the state-based rows are flat — no module-grain cell exceeds
|r| = 0.35.

### Caveat recorded in `metrics_lib.R`

`T_inter_aligned` is **signed**, and its sign flips if a module's two teams
are relabelled. The orientation is pinned by the `TEAM1` table, so values are
internally consistent, but the sign is only meaningful relative to that
convention. `|T_inter_aligned|` and `T_inter_blocks` are convention-free and
should be preferred unless the orientation is being asserted on purpose. This
matters for reading the strongly negative "inter-module alignment" column in
matrices 2-5.


---

## Figure 1 rebuilt on RACIPE (2026-09-24)

**Output folder.** The revised Figure 1 writes to **`Figures/figure1_new/`**,
separate from `Figures/panels/`, which keeps the earlier FIGURE_PLAN.md
layout (Ising Fig 1, plus Figs 2-4). The two are different designs of the same
figure and should not be mixed: `Figures/panels/fig1b_module_properties.png`
and `Figures/figure1_new/fig1C_module_properties.png` are the old and new
versions of the same content. Nothing in `Figures/panels/` was changed.

Directive change: **Figure 1 leads with RACIPE** — fewer metrics and less
ambiguity than the Ising side — and the argument runs PC metrics first, then
P_s. This reverses the 2026-09-22 "Ising only, RACIPE for PCA only"
convention *for these panels*; downstream P_s panels stay Ising, and any
figure mixing the two must say so on its face.

Source outline: `Figures/correlation_heatmaps/figure_outline_new.md`.

### Getting the data off the cluster

The RACIPE solutions live on Explorer at
`/scratch/a.hari/InteractingTeams/<order>_EMT_MechSens/RACIPE/` — **1.9 GB for
the doubles alone**, so nothing bulky is transferred. `racipe_z.sh` (copied to
the cluster) does the work remotely in one awk pass per network:

1. pools the three replicates `<net>_{A,B,C}_solution.dat`;
2. z-scores each gene column across the pooled set;
3. emits a systematic subsample (default 3000 rows).

Only the subsample comes back, into `RACIPE_zscores/<net>_z.csv`. Total
transferred: **1.8 MB for seven networks**, against ~150 MB of source.

Two things about the cluster files that the local `6_module_racipe/` copies
lack: they carry a **header row with gene names**, and a `basin` column. The
layout is `parIndex nStates basin <genes...>`, one row per state. Replicates
are suffixed `_A/_B/_C`, and the **self-pairs `Comb_<M>_<M>` are present**, so
single-module RACIPE needs no new runs.

### Panels built

`Figures/fig1_racipe.R` -> `Figures/figure1_new/`:

| file | content |
|---|---|
| `fig1AB_racipe_states.png` | A: RS (T_s 0.92), B: PS (T_s 0.09) |
| `figS_racipe_all_modules.png` | all six modules, weakest to strongest team |
| `fig1E_racipe_PS_RS.png` | PS + RS coupled — the same pair as A/B |

Nodes are ordered by team from the `allNodes.nodes` `NewProgram` split, with a
dashed line at the team boundary and a solid line at module boundaries. Rows
are ordered by the team-1 minus team-2 contrast rather than clustered, so the
bimodality is visible without imposing a clustering choice. Colour is squished
at +/-2.5 z.

**The panels do the job A/B were meant to do.** RS shows two crisp mirrored
team blocks with a sharp switch; PS is mottled with no block structure. The
six-module supplementary reads as a clean progression with team strength:
PS mottled, C and M intermediate, E / AS / RS sharply bimodal.

### Panels C-G (2026-09-24)

`Figures/fig1_panels_CDEFG.R` -> `Figures/figure1_new/`.

| panel | file | note |
|---|---|---|
| C | `fig1C_module_properties.png` | team strength, within-module density, cross-module impurity, as lollipops ordered by team strength |
| D | `fig1D_pca_correlations.png` | r with 95% CI against the three PCA metrics |
| E | `fig1E_influence_PS_RS.png` | influence matrix of PS+RS, **right third left blank for manual annotation** |
| F | `fig1F_intra_Ts_by_nmodules.png` | intra-module T_s, violins per module dodged by module count |
| G | `fig1G_pc_vs_nmodules.png` | the three PC metrics against module count |

**C** uses *cross-module* impurity, because within-module impurity is 0 for
five of the six modules and has no variance to plot.

**E** is the bare matrix with module boundaries (solid), team splits (dashed)
and per-module T_s labels; `plot_spacer()` reserves the right third.

**G** z-scores each PC metric so three different units share one axis — never
a second y-axis. The trend is monotone and clean:

| modules | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| PC1 variance (z) | +2.15 | +0.48 | -0.25 | -0.66 | -0.88 |
| PCs to 90% (z) | -1.50 | -0.67 | +0.03 | +0.71 | +1.30 |
| PCA entropy (z) | -1.50 | -0.26 | +0.20 | +0.40 | +0.49 |

Adding modules monotonically flattens the RACIPE spectrum — the cleanest
"diversity increases with coupling" statement in the project, and it uses no
phenotypic score at all.

### A third Simpson's paradox, found while building D

`T_network_intra_sum` is a **sum over modules**, so it is +0.837 correlated
with module count. Its pooled correlation with PC1 variance is **-0.571**
while every within-module-count value is **positive** (+0.34, +0.31, +0.15,
+0.13). The mean is deconfounded (r = -0.03 with module count) and pooled
+0.225, consistent with the within-count signs. **Panel D now uses the mean.**

This is the same defect already found in `T_inter_sum` (inter-module) — any
metric summed over modules or module pairs inherits module count. Worth a
sweep: `T_network_intra_sum` and `T_inter_sum` both appear as inputs in the
correlation matrices and in `loo_results.md`.

To keep this visible rather than merely caveated, **panel D plots two points
per cell**: filled = pooled r with CI, hollow = mean of the within-count
correlations. Where they straddle zero the pooled value is a module-count
artefact. Three cells currently flip, all for inter-module T_s (against PC1
variance, PCs to 90%, and PCA entropy).

Other confounded predictors, flagged in the caption rather than removed:
network size (r = +0.94 with module count) and density (-0.80). Both keep
their sign within counts, so they are amplified rather than inverted.

### Still to build

- Supplementary row 1: networks and interaction matrices for A and B.
- Final assembly into the outline's row layout (A / B / C+D / E | F+G).


## Open decisions

1. **How much of Figure 1B survives?** Resolved for 1e/1f/1j — all three are
   built. 1e (frustration, r = -0.816) earns a main-text slot; 1f/1j go to
   supplementary as a methods caveat until `F_hy` is redefined. Still open:
   **is 1g (`P_mean`) worth a cluster run?** It is the only genuinely new
   simulation left in Figure 1.
2. **Figure 5 before or after Figure 4's verdict?** Recommendation: after. If
   Figure 4 identifies no predictive metric, redirect the effort into expanding
   the corpus — giving two-team splits to some of the 8 excluded single-team
   modules would take n from 6 toward 14 and make the property regressions in
   Q2/Q4 meaningful for the first time.
3. **Leave out module *pairs* as well as singles?** 15 pairs + 6 singles = 21
   removal effects instead of 6, regressed against interaction properties
   rather than module properties. `metrics_pair_level.csv` already holds the
   predictors. This is the cheapest route to making Figure 4 positive.
4. **RS vs PS or AS vs PS** as the strong/weak contrast in 1a.
8. **Adopt the team-aware inter-module metric, and mean- not sum-normalise?**
   Both changes are well-supported (above) and together reverse Figure 4b's
   sign. The cost is that the published `inter_module_sum_T_s` rows no longer
   correspond. Recommend adopting `T_inter_blocks_mean` and keeping
   `T_inter_sum` in supplementary for continuity.
7. **How is the matrix-weighted P_s framed in the text?** It is a
   team-structure alignment measure (`|w^T M S|`), a team-resolved sibling of
   frustration — not an independent landscape readout. Options: present
   Figure 4 as an *alignment* result (defensible as it stands), add the
   synthetic-committed baseline as the comparison, or lead with the
   state-based score. This changes Figure 4's claim, not its data.
6. **Which state-based P_s definition is canonical** — the stated one (used
   here) or the published `_stateps` lineage with its per-network `1/N`
   factor? They give different correlations, not just different scales. The
   `dominant` rows reconcile exactly under `/N`; `top30`/`top60` never do.
   Worth settling before any state-based number goes in the paper.
5. **`dominant` vs `top60`** for the Fig 2b/2c heatmaps — `dominant` is a
   coin toss in 5 of the 56 networks, all two-module, which is exactly the
   subset Figure 2 uses.
