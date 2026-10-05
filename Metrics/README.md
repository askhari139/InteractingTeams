# Metrics from the formula sheet

R implementation of the metrics in
`Tkl_component_scatterplots/StateBased_newnorm/lmax10/outline.md` (the formula
image), built on `funcsKishore`.

| file | what |
|---|---|
| `FORMULAE.md` | **every formula as implemented**, with the departures from the formula image called out — read this to audit the maths |
| `metrics_lib.R` | the metric functions |
| `compute_metrics.R` | driver; writes the two CSVs below |
| `validate_metrics.R` | checks topology metrics, `T_network`, `T_module_intra` |
| `validate_ps_inter.R` | checks inter-module `T_A->B`; prints P_s comparison |
| `validate_state_ps.R` | checks the state-based phenotypic score |
| `validate_infl_ps.R` | checks the influence-weighted phenotypic score |
| `validate_network_ps.R` | checks the network-level phenotypic score |
| `metrics_network_level.csv` | network x l_max |
| `metrics_module_level.csv` | network x module x l_max |
| `metrics_pair_level.csv` | network x (module A -> module B) x l_max |

Run from the repo root:

```sh
Rscript Metrics/compute_metrics.R          # default l_max = 1,5,7,10
Rscript Metrics/compute_metrics.R 1,10     # or pick your own
Rscript Metrics/validate_metrics.R
```

## Inputs

- `<net>.topo`, `<net>.teams` — network topology and its 2-team partition
  (the `.teams` file labels only *core* nodes; peripherals are pruned before
  team detection).
- `Single_<M>.topo`, `Single_<M>_nodes.txt` — per-module topology and the node
  order used by the state strings.
- `allNodes.nodes` — **the authoritative module-level team definition.**
  `Program` is the module, `NewProgram` splits it into its two teams. Do *not*
  use the `Single_<M>.teams` files for module teams: they disagree with this
  for PS, and only the `NewProgram` split reproduces the published numbers.
- `<net>_finFlagFreq.csv` — Boolean steady states with basin sizes (`Avg0`);
  in `Two_module_nets/` for 2-module networks, in the per-order folders
  otherwise. The `states` string is ordered as `<net>_nodes.txt`.

## Verified against the existing `correlations_*` outputs

Every check below is exact (< 1e-9) against
`Combined_data/combined_module_level_data.csv` and the `*_panel_*.csv` files.

| metric | result |
|---|---|
| `density` | 56 / 56 networks |
| `impurity` | 56 / 56 |
| `intra_density` | 56 / 56 |
| `inter_density` | 56 / 56 |
| `T_network` | 224 / 224 (56 nets x 4 l_max) |
| `T_A->B` (inter-module) | **1800 / 1800** |
| `T_module_intra` | **720 / 720** (all modules including PS) |
| state-based `P_s_dominant` | 170 / 180 (the 10 exceptions are basin ties, below) |
| influence-weighted `P_s_dominant` (module level) | **720 / 720** |
| network-level `P_s_dominant` (via published correlations) | **80 / 80** |

## Where the implementation departs from the formula image

Four places (the fourth, the `+/-1` coding of the state-based score, is in
`FORMULAE.md` §6: the implemented `P_s` is exactly `2 x (E(T_1) - E(T_2))`,
so its range is `[-2, 2]` rather than the sheet's `[-1, 1]`. Correlations are
unaffected; absolute thresholds are not). In each case the code follows what reproduces the published
numbers, and exposes the literal reading as an option.

1. **`T_kl` takes the absolute value after the block mean, not before.** The
   sheet writes `T_kl = (1/|E_kl|) * sum |I_ij|`. The pipeline computes
   `|mean(I_ij)|`, so opposing edges cancel. This also matches
   `funcsKishore::getGsVec()`. `Tkl(..., abs_first = TRUE)` and
   `T_inter(..., abs_first = TRUE)` give the sheet's version. The difference is
   large: for `Comb_AS_C` at l_max 10, `T_AS->C` is 0.0282 (implemented) vs
   0.0508 (literal).

2. **`|E_kl|` is all node *pairs*, not edges.** Read literally as edges, every
   `T_kl` would be exactly 1 at l_max = 1, since the influence matrix is just
   the signed adjacency there. The denominator is `|team_k| * |team_l|`.

3. **`Density` uses `N(N-1)`, not `N^2`,** over the core subgraph with
   self-loops dropped. This is what reproduces all 56 networks; `N^2` does not.

Two further conventions worth recording, both needed to reproduce the outputs:

- **The inter-module metrics need the *unpruned* influence matrix.** In
  `Comb_AS_C` the single AS->C edge (`Casp3 -| N_bcatenin`) lands on a
  peripheral node, so `funcsKishore::InfluenceMatrix()`'s reduced matrix zeroes
  it out. `influence_full()` reproduces `InfluenceMatrix()` exactly on the core
  block but keeps peripherals.
- **`T_module_intra` uses the module's *own* teams on the *network's*
  influence matrix.** That is why it is network-independent at l_max = 1 (no
  path leaves the module) but network-dependent at l_max >= 5 (paths run
  through the other modules) — exactly the pattern in the published data.
  The `NewProgram` labels cover every module node, peripherals included, so no
  team-assignment heuristic is needed.

## Module teams, and the P_s sign convention

Only six modules have a two-team split in `allNodes.nodes`:

| module | Program | team 1 (positive pole) | team 2 |
|---|---|---|---|
| AS | Apoptotic_Switch | Apoptotic_Pro (13) | Apoptotic_Anti (4) |
| C | CIP | CIP_E (6) | CIP_M (8) |
| E | EMT | EMT_Mes (9) | EMT_Epi (3) |
| M | Migration | Migration_Pro (5) | Migration_Anti (1) |
| PS | Phase_Switch | Phase_Switch (12) | Phase_acyclic (3) |
| RS | Restriction_Switch | Restriction_Cyclin (4) | Restriction_arrest (3) |

`A`, `CC`, `GBM`, `GP`, `Gc`, `Gm`, `OL` and `Tbs` have a single `NewProgram`,
so no team strength or phenotypic score is defined for them. `compute_metrics.R`
therefore skips any network containing one — which is exactly why the published
`correlations_*` analysis covers 56 networks over 6 modules and excludes CC.

`P_s = mean(S over team1) - mean(S over team2)`, so the `TEAM1` table in
`metrics_lib.R` fixes the sign of every phenotypic score. The orientations above
reproduce the published signs (AS 30/30; C, M, PS, RS 28/30). **EMT is oriented
mesenchymal-positive**; taking `EMT_Epi` first flips the sign in 28 of 30
networks. Change `TEAM1` if you want a different convention.

## Known gaps

1. **`P_s_dominant` is ambiguous in 5 networks.** 170/180 module-network pairs
   reproduce the published value exactly. All 10 exceptions are exact sign
   flips in `Comb_C_E`, `Comb_C_M`, `Comb_E_RS`, `Comb_M_PS` and `Comb_PS_RS`,
   and in every one the top two states are near-tied mirror attractors
   (0.27662 vs 0.27638; 0.45273 vs 0.45234; 0.20016 vs 0.19884; 0.11036 vs
   0.11011; 0.16535 vs 0.16385). Which one is "dominant" turns on the 3rd-4th
   decimal of the basin size, so it is not robust to simulation noise. Prefer
   the `top30`/`top60` aggregates over `dominant` for these networks.

2. **The published `_stateps` columns carry an extra `1/N_network` factor.**
   The current definition takes a mean difference and needs no further
   normalisation, and `dominant` reproduces as `published * N_network`. The
   `top30`/`top60` columns do not reproduce under any constant factor, so the
   old aggregates were computed differently; the implementation here follows
   the stated definition (weighted *sum* over the top 30%/60% of the state
   space) rather than the old columns.

3. **The basin frequencies have drifted from the published run.** The
   influence-weighted `P_s_dominant` now reproduces 720/720, but `top30` and
   `top60` split perfectly along one line (the same holds for the state-based
   variant: 2/2 reproduced where weights are irrelevant, 0/358 where they
   matter, with the published value inside the achievable range in all 358):

   | selected states | reproduced exactly |
   |---|---|
   | all share one \|dF\| (weights irrelevant) | 223/223 top30, 59/59 top60 |
   | \|dF\| varies (weights matter) | 0/497 top30, 0/661 top60 |

   In every non-reproducing case the published value still lies inside
   `[min|dF|, max|dF|]` of the same selected states (497/497 and 661/661), so
   it is a weighted average of identical per-state scores with *different*
   weights. Median relative error 1.0% (top30) and 0.34% (top60). The `Avg0`
   values in the current `finFlagFreq.csv` files are therefore not the ones the
   published columns were computed from -- expected Monte-Carlo scatter between
   Boolean runs. This is also why 10 state-based `dominant` values flip sign:
   all are near-tied mirror attractors.

4. **The two P_s variants use different aggregations.** `P_s_state_*` uses a
   signed, unnormalised weighted *sum* (the stated definition).
   `P_s_infl_*` uses `|dF|` with renormalised weights (a weighted *mean*),
   because that is what reproduces the published column -- a weighted sum would
   be ~0.3x smaller and would fall outside `[min|dF|, max|dF|]`, which it never
   does. Pass `absolute=` / `normalise=` to `aggregate_ps()` to align them.

## Caution when using funcsKishore

`findClusterTeams()` **overwrites `<net>.teams`** (`writeLines(l, paste0(net,
".teams"))`). Do not call it on a network whose `.teams` file you want to keep.
`InfluenceMatrix()` writes into `Influence/` and `setwd()`s during the call.

`getGsVec()` does **not** produce the PS 12/3 split. With `method = "Cluster"`
it reads an existing `.teams` file (so it returns the 6/9 partition, Gs =
0.1779) and only calls `findClusterTeams()` when that file is absent;
`method = "Brute"` gives Gs = 0.2931. Neither equals the published PS intra
value of 31/144 = 0.2153. Only the `allNodes.nodes` `NewProgram` split does.
Note also that `getGsVec()`'s `Gs` is a *different* quantity from `T_network`:
it builds cross-team sub-topologies and computes the influence matrix on each,
whereas `T_network` uses the whole-network influence matrix.
