# Leave-one-module-out effect on correlations

Answers Q1 and Q2 of
`Tkl_component_scatterplots/StateBased_newnorm/lmax10/analysis.md`, following
the method in its line 10. Q3 and Q4 (clustering inconsistency) are **not**
addressed here — they need a clustering analysis rather than a correlation one.

Scripts: `module_properties.R`, `loo_analysis.R`, `loo_plots.R`.
Outputs: `module_properties.csv`, `loo_correlations.csv`, `loo_null.csv`,
`loo_impact_vs_properties.csv`, `figures/fig1..3*.png`.

## Method

- Network-level P_s is the **mean of its module P_s values** — the aggregation
  that reproduces the published `*_module` correlation rows (80/80 on
  `dominant`), so it is the project's own convention rather than a new choice.
- `r` = Pearson correlation of that P_s against `T_network`, over networks.
- Baseline: all 56 networks. Leave-one-out: for module X, only the 26 networks
  that do not contain X — the same construction the `_noPS` folders use for PS.
  Every module leaves exactly 26 networks, so the removals are balanced.
- Impact = `r_all − r_without_X`.
- Four P_s measures: state-based and influence-weighted, `top30` and `top60`,
  at each of l_max 1, 5, 7, 10.
- **Null model:** `r` recomputed on 2000 random 26-of-56 subsets per measure
  per l_max. Without this the impacts are unreadable — `r` on 26 points is
  noisy, and the question "which module matters most" is meaningless until the
  subsampling spread is known.

## Baseline correlations (all 56 networks)

Range across l_max = 1, 5, 7, 10:

| predictor | outcome family | r |
|---|---|---|
| network T_s | network P_s, influence-weighted | +0.930 to +0.954 |
| network T_s | module P_s, influence-weighted | +0.900 to +0.936 |
| network T_s | PCA metrics | −0.850 to +0.823 |
| network T_s | network P_s, state-based | −0.199 to +0.367 |
| network T_s | module P_s, state-based | −0.030 to +0.229 |
| inter-module T_s sum | PCA metrics | −0.384 to +0.634 |
| inter-module T_s sum | module P_s, influence-weighted | −0.567 to −0.080 |
| inter-module T_s sum | network P_s, influence-weighted | −0.561 to −0.211 |
| inter-module T_s sum | module P_s, state-based | −0.136 to +0.029 |
| inter-module T_s sum | network P_s, state-based | −0.051 to +0.141 |

Four things stand out before any module is removed:

1. **Network-level P_s is the strongest and cleanest readout** (r = 0.93–0.95
   against network T_s) — slightly better than the module-level
   influence-weighted P_s, and with a much tighter null (below).
2. **The state-based P_s barely correlates with team strength at all**
   (|r| <= 0.23 against network T_s, <= 0.14 against the inter-module sum).
   The two P_s variants are not interchangeable readouts.
3. **PCA metrics track network T_s strongly and with the expected signs:**
   PC1 variance +0.80 to +0.61 (falling with l_max), PCs-to-90% −0.76 to −0.40,
   PCA entropy −0.85 to −0.58. Stronger team strength means a more concentrated,
   lower-entropy landscape.
4. **Inter-module T_s sum correlates *negatively* with every P_s measure**, and
   the correlation weakens sharply with l_max (e.g. network-level P_s: −0.56 at
   l_max 1, −0.21 at l_max 10). More inter-module coupling goes with a weaker
   phenotypic score — the opposite sign to network T_s.

## The four P_s variants

Both levels now carry both formulae — the full 2x2:

| | influence-weighted | state-based |
|---|---|---|
| **module-level** | `P_s_infl_*`, averaged over the network's modules | `P_s_state_*`, same aggregation |
| **network-level** | `P_s_net_infl_*`, the network's own teams on the core matrix, \|dF\| | `P_s_net_state_*`, signed team-mean difference |

The module-level aggregation is the mean over the network's modules — the
convention that reproduces the published `*_module` correlation rows (80/80 on
`dominant`). The network-level influence-weighted score is the quantity behind
the published `*_network` rows (also 80/80).

**Influence-weighted beats state-based at both levels, and network beats module
within the influence-weighted pair.** Against network T_s the ordering is
network-infl (0.93-0.95) > module-infl (0.90-0.94) >> network-state (-0.20 to
+0.37) > module-state (-0.03 to +0.23). The influence-weighted scores are also
far more stable: their nulls are ~0.045 against ~0.33-0.45 for the state-based
ones, so the state-based correlations are essentially unusable as a readout at
this sample size.

One deliberate asymmetry: the network-level influence-weighted score is `|dF|`
because that is what reproduces the published coefficients, while every other
variant is signed.

## Q1 — which module has the highest effect on correlations?

Mean |Δr| over the four l_max, per predictor and outcome family. The
highest-impact module is **bold**:

| predictor | family | ranking |
|---|---|---|
| network T_s | module P_s, infl | **PS 0.047**, M 0.035, AS 0.021, C 0.019, RS 0.013, E 0.007 |
| network T_s | module P_s, state | **E 0.178**, M 0.175, RS 0.134, C 0.118, AS 0.096, PS 0.068 |
| network T_s | network P_s, infl | **PS 0.017**, AS 0.016, E 0.016, M 0.014, RS 0.013, C 0.013 |
| network T_s | network P_s, state | **C 0.288**, RS 0.123, E 0.108, M 0.073, PS 0.072, AS 0.050 |
| network T_s | PCA metrics | **RS 0.148**, AS 0.142, M 0.103, PS 0.085, C 0.073, E 0.038 |
| inter-module T_s sum | module P_s, infl | **M 0.249**, AS 0.159, RS 0.158, PS 0.100, C 0.095, E 0.050 |
| inter-module T_s sum | module P_s, state | **E 0.122**, RS 0.107, C 0.093, PS 0.059, M 0.049, AS 0.024 |
| inter-module T_s sum | network P_s, infl | **M 0.217**, AS 0.158, RS 0.136, PS 0.059, C 0.043, E 0.042 |
| inter-module T_s sum | network P_s, state | **C 0.180**, E 0.083, RS 0.075, PS 0.070, M 0.068, AS 0.040 |
| inter-module T_s sum | PCA metrics | **AS 0.141**, M 0.113, C 0.097, E 0.092, RS 0.088, PS 0.033 |

**There is no single answer. All six modules lead at least one cell** — PS, E,
C, RS, M and AS each top some correlation. That inconsistency is the result.

Against the null, 38 of 576 removals exceed the random-subset 95th percentile
where chance gives ~29. Per family:

| predictor | family | null p95 | max \|Δr\| | beyond null |
|---|---|---|---|---|
| network T_s | network P_s, infl | 0.045 | 0.034 | 1 / 48 |
| network T_s | module P_s, infl | 0.046 | 0.074 | 7 / 48 |
| network T_s | module P_s, state | 0.334 | 0.392 | 3 / 48 |
| network T_s | network P_s, state | 0.453 | 0.393 | 0 / 48 |
| network T_s | PCA metrics | 0.205 | 0.278 | 8 / 96 |
| inter-module T_s sum | module P_s, infl | 0.210 | 0.385 | 8 / 48 |
| inter-module T_s sum | network P_s, infl | 0.204 | 0.328 | 6 / 48 |
| inter-module T_s sum | module P_s, state | 0.275 | 0.187 | 0 / 48 |
| inter-module T_s sum | network P_s, state | 0.287 | 0.210 | 0 / 48 |
| inter-module T_s sum | PCA metrics | 0.203 | 0.244 | 5 / 96 |

Reading this table:

- **Against network T_s, nothing moves.** The two strong readouts have nulls of
  ~0.045 and impacts that barely clear them (network-level P_s: 1 of 48). The
  correlation is a property of the corpus, not of any one module.
- **The state-based P_s has a null of 0.27–0.34** — a third of the correlation
  range — so its large-looking impacts (up to 0.39) carry no information. Its
  baseline r is ~0, which is why any resampling swings it.
- **The inter-module T_s sum is where module identity matters most**: M's
  removal shifts that correlation by up to 0.385 (module P_s) and 0.328
  (network P_s), and 8 of 48 and 5 of 48 removals clear the null. If any part
  of this analysis is worth pursuing, it is M against the inter-module sum.
- **PS leads against network T_s but with tiny absolute effects** (0.047 and
  0.017) — and it is *last* for the inter-module sum with PCA metrics. Its
  special treatment in the `_noPS` folders is not reflected in a distinctive
  effect on these correlations.

## Q2 — module properties, and which predicts the effect best

`module_properties.csv`, at l_max 10:

| module | team strength | density within | impurity within | intra-team density | inter-team density | density across | impurity across |
|---|---|---|---|---|---|---|---|
| AS | 0.766 | 0.208 | 0.000 | 0.167 | 0.271 | 0.0103 | 0.368 |
| C | 0.497 | 0.194 | 0.000 | 0.281 | 0.125 | 0.0144 | 0.957 |
| E | 0.715 | 0.344 | 0.000 | 0.293 | 0.438 | 0.0071 | 0.600 |
| M | 0.576 | 0.350 | 0.000 | 0.417 | 0.250 | 0.0244 | 0.842 |
| PS | 0.086 | 0.300 | 0.143 | 0.304 | 0.292 | 0.0167 | 0.357 |
| RS | 0.917 | 0.405 | 0.000 | 0.333 | 0.333 | 0.0279 | 0.360 |

Two things to flag about the properties themselves:

- **`impurity_within` is 0 for five of six modules** (only PS is non-zero, at
  0.143). It has essentially no variance and cannot serve as a predictor. PS
  being the sole impure module — and by far the weakest team, 0.086 against
  0.50–0.92 for the rest — is the cleanest structural fact in this table.
- **Cross-module impurity is high and variable** (0.36–0.96). C and M connect to
  other modules almost entirely through team-misaligned edges.

Spearman correlation of mean |Δr| against each property, averaged over all
measures and l_max, per predictor (`loo_impact_vs_properties.csv`):

| property | vs network T_s | vs inter-module T_s sum |
|---|---|---|
| density across modules | +0.16 | +0.09 |
| impurity across modules | −0.15 | +0.14 |
| density within | +0.10 | +0.12 |
| team strength | +0.09 | +0.17 |
| inter-team density (within) | +0.06 | −0.06 |
| intra-team density (within) | +0.05 | +0.05 |

**No property predicts the effect, under either predictor.** Every mean |ρ| is
below 0.2, and `impurity_across` even flips sign between the two predictors.
Individual panels in `properties_*.png` reach ρ = +0.83 (PCA metrics vs density
across) and ρ = +0.77 (state-based P_s vs impurity across), but with n = 6 those
are not significant and the same property lands near zero in adjacent panels.

The honest answer to "which property measures the effect best" is that **this
design cannot tell**: six modules is too few, and (Q1) most of the impacts are
resampling noise to begin with.

## What would make this answerable

1. **Use the influence-weighted P_s as the readout.** Its baseline r ≈ 0.93 and
   its null spread is ~5× tighter than the state-based one, so real effects have
   somewhere to show up.
2. **Get more modules.** The 8 single-team modules (A, CC, GBM, GP, Gc, Gm, OL,
   Tbs) are excluded because `allNodes.nodes` gives them one `NewProgram`. If
   any can be given a two-team split, n goes from 6 toward 14 and the property
   regressions become meaningful.
3. **Leave out module *pairs* as well as singles.** 15 pairs plus 6 singles
   gives 21 removal effects to regress against interaction properties, and the
   pair-level table (`metrics_pair_level.csv`) already holds the predictors.
4. Treat `r` differences against a resampling null throughout — `loo_null.csv`
   has the machinery.

## Figures

All in `Metrics/figures/`, styled with `funcsKishore::theme_Publication`.

| file | content |
|---|---|
| `violin_<predictor>_<family>.png` (10) | Distribution of correlations. Grey violin: r over 2000 random 26-of-56 subsets. Blue points: the six leave-one-module-out values, direct-labelled. Orange diamond: the original r over all 56 networks. x = l_max. |
| `impact_<predictor>_<group>.png` (6) | Δr by module, one point per l_max. Point fill is the exact correlation the removal produced, on a fixed −1…+1 diverging scale with an outline so near-zero r stays visible. Grey band: central 90% of the null. |
| `properties_<predictor>_<group>.png` (6) | Δr against each module property, coloured by module, faceted family × property, Spearman ρ annotated. |

`<predictor>` is `T_network` or `T_inter_sum`; `<group>` is `psM` (module-level
P_s), `psN` (network-level P_s) or `pc` (PCA metrics); `<family>` splits the
P_s groups further into influence-weighted and state-based.

Families are plotted separately because their correlation scales differ by an
order of magnitude — pooling the state-based P_s (r ~ 0) with the
influence-weighted one (r ~ 0.93) in a shared axis would hide exactly the
difference that matters.

Module colours are the 6-of-8 subset of the reference categorical palette that
clears the all-pairs colourblind (ΔE 6.9) and normal-vision (ΔE 15.6) floors;
since ΔE 6.9 sits in the band that requires secondary encoding, every module
point also carries a direct text label.
