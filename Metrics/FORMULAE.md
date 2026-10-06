# Formulae as actually implemented

Every formula below is transcribed from `Metrics/metrics_lib.R`, function by
function. Where the implementation differs from the formula image
(`Tkl_component_scatterplots/StateBased_newnorm/lmax10/outline.md`) the
difference is called out in a **Departure** note, with the reason. There are
four: the three in `T_kl` / `|E_kl|` / `Density` (§3, §5), and the `+/-1`
coding that doubles the state-based phenotypic score (§6).

The right-hand column of each heading gives how the quantity was checked
against the existing `correlations_*` outputs.

---

## 0. Notation and index conventions

| symbol | meaning |
|---|---|
| `V` | all nodes of the network, in the order given by `<net>_nodes.txt` |
| `N` |  number of nodes in the **whole network** |
| `A` | signed adjacency: `A[i,j] = +1` if `i -> j` activates, `-1` if it inhibits, `0` otherwise |
| `I` | influence matrix (§1) |
| `M` | the node set of one module, from `Single_<M>_nodes.txt` |
| `T_1, T_2` | the two teams of a partition |
| `S` | a Boolean steady state in ±1 coding |
| `f` | basin size of a state (the `Avg0` column of `finFlagFreq.csv`) |

**Direction convention, used everywhere:** the first index is the *source*,
the second the *target*. `A[i,j]` and `I[i,j]` mean *the effect of i on j*.
This follows `funcsKishore::TopoToIntMat`, which fills `intmat[source, target]`.
Getting this backwards changes the inter-module and phenotypic-score results
(the matrices are not symmetric), so it is worth checking first.

`.topo` files code `Type = 1` as activating and `Type = 2` as inhibiting; both
`read_topo()` and `TopoToIntMat()` map `2 -> -1`.

---

## 1. Influence matrix — `influence_full()`

For a path length cutoff `l_max`:

```
I = (1 / l_max) * sum_{l=1}^{l_max} R_l

           A^l [i,j]
R_l[i,j] = -----------      when |A|^l [i,j] != 0
           |A|^l [i,j]

R_l[i,j] = A^l [i,j]        when |A|^l [i,j] == 0   (the 0/0 case)
```

`A^l` is the ordinary matrix power, so `A^l[i,j]` is the signed sum over all
length-`l` paths from `i` to `j` and `|A|^l[i,j]` is the number of such paths.
Each `R_l` is therefore a per-path-length *consistency* in `[-1, 1]`, and `I`
averages that over `l = 1 .. l_max`.

Two consequences worth keeping in mind while auditing:

- **At `l_max = 1`, `I = A` exactly** (`R_1 = A/|A| = sign(A)`, and entries with
  no edge fall into the 0/0 branch and become `0`). Every `l_max = 1` number
  below is therefore a pure edge count, which is what made the conventions
  checkable by hand.
- `|I[i,j]|` grows with `l_max` in general, so every quantity built from `I` is
  `l_max` dependent.

**Departure — peripherals are kept.** `funcsKishore::InfluenceMatrix()`
iteratively prunes *peripheral* nodes (those that are not both a source and a
target) and returns the reduced matrix. `influence_full()` runs the identical
recursion but keeps all `N` rows and columns. On the core sub-block the two
agree to 0 exactly (checked). Pruning cannot be used here because
inter-module edges can land on peripheral nodes — in `Comb_AS_C` the *only*
AS->C edge is `Casp3 -| N_bcatenin` and `N_bcatenin` is peripheral, so the
reduced matrix would report `T_AS->C = 0` instead of `1/238`.

---

## 2. Teams

### 2.1 Network teams

Read verbatim from `<net>.teams`: one comma-separated team per line, two lines.
These files label **only core nodes** — peripherals are pruned before team
detection, so `unlist(teams)` is the core node set. `Comb_AS_C.teams` covers 25
of 31 nodes.

### 2.2 Module teams — `module_teams()`

Taken from `allNodes.nodes`, **not** from the `Single_<M>.teams` files:
`Program` names the module and `NewProgram` splits it in two.

| module | Program | team 1 | team 2 |
|---|---|---|---|
| AS | Apoptotic_Switch | Apoptotic_Pro (13) | Apoptotic_Anti (4) |
| C | CIP | CIP_E (6) | CIP_M (8) |
| E | EMT | EMT_Mes (9) | EMT_Epi (3) |
| M | Migration | Migration_Pro (5) | Migration_Anti (1) |
| PS | Phase_Switch | Phase_Switch (12) | Phase_acyclic (3) |
| RS | Restriction_Switch | Restriction_Cyclin (4) | Restriction_arrest (3) |

`A`, `CC`, `GBM`, `GP`, `Gc`, `Gm`, `OL`, `Tbs` have a single `NewProgram` and
so no two-team structure; networks containing them are skipped entirely.

**Why not `Single_<M>.teams`:** those files disagree for PS — they give 6/9
(`Cdc20, Cdh1, Emi1, UbcH10, Wee1, pAPC` vs 9) where `NewProgram` gives 12/3.
Only the 12/3 split reproduces the published PS intra-module value of `31/144`;
an exhaustive scan of all 16384 partitions of PS's 15 nodes found it to be the
unique partition that does. The two sources agree for the other five modules.

**Team order fixes every P_s sign** (§6 and §7 are differences `T_1 - T_2`).
The order in the table above reproduces the published signs: 30/30 for AS,
28/30 for C, M, PS, RS. **EMT is oriented mesenchymal-positive** — putting
`EMT_Epi` first flips the sign in 28 of 30 networks. This is the one place a
convention was chosen rather than derived, so it is the first thing to check if
a P_s sign looks wrong. It is the `TEAM1` constant in `metrics_lib.R`.

---

## 3. Team strength

### 3.1 `T_kl` — `Tkl()`

For teams `T_k`, `T_l`:

```
             |    1                               |
T_kl =       | -------- *  sum       sum   I[i,j] |
             | |T_k||T_l|   i in T_k   j in T_l   |
```

i.e. **the absolute value of the mean of `I` over the block**.

**Departure — two of them, both needed to reproduce the outputs.**

1. *The absolute value is taken after the mean, not before.* The sheet writes
   `T_kl = (1/|E_kl|) sum |I_ij|`. Under the sheet's reading opposing edges
   would not cancel; under the implemented reading they do. This matches
   `funcsKishore::getGsVec()`, which computes signed block sums and takes
   `abs()` afterwards. The gap is not small: for `Comb_AS_C` at `l_max = 10`,
   `T_AS->C` is `0.0282` as implemented against `0.0508` read literally.
   `Tkl(..., abs_first = TRUE)` gives the literal version.
2. *`|E_kl|` is the number of node **pairs**, `|T_k| * |T_l|`, not the number of
   edges.* Read as edges, every `T_kl` would be exactly `1` at `l_max = 1`,
   since `I = A` there and `|A[i,j]| = 1` on every edge. The published
   `l_max = 1` values are not 1, so the denominator must be all pairs.

### 3.2 Network team strength — `T_network()`  — **224/224 exact**

Over the network's own two teams:

```
T_network = (T_00 + T_01 + T_10 + T_11) / 4
```

The mean is over the four blocks, equally weighted — not weighted by block
size. Teams are first intersected with `rownames(I)`.

### 3.3 Intra-module team strength — `T_module_intra()` — **720/720 exact**

Same four-block average, but the teams are **the module's own two teams**
(§2.2) while the matrix is **the whole network's `I`**:

```
T_M^mod = mean over the 4 blocks of  T_kl  ,  with  T_k, T_l  the teams of M
```

This mixture is the part most worth scrutinising, and it is also what the
published data demands. At `l_max = 1` no path can leave the module, so the
within-module block of the network's `I` equals the module's own `A` and
`T_M^mod` is network-independent — and indeed the published values at
`l_max = 1` take exactly one distinct value per module. At `l_max >= 5` paths
run out through the other modules and back, so `T_M^mod` becomes
network-dependent — and the published values then take 2-21 distinct values per
module. Using the module's *isolated* topology instead would keep it constant
at all `l_max`, contradicting the data.

### 3.4 Inter-module team strength — `T_inter()` — **1800/1800 exact**

For modules `A`, `B` (all their nodes, not their teams):

```
              |    1                            |
T_{A->B} =    | ------- *  sum      sum   I[i,j] |
              | |A| |B|   i in A   j in B        |
```

Absolute value after the mean, as in §3.1; the denominator is `|A| * |B|`, all
ordered pairs. Computed on the unpruned `I` (§1).

Worked check at `l_max = 1`, where `I = A`: `Comb_AS_C` has exactly one
cross-module edge, `Casp3 -| N_bcatenin`, with `|AS| = 17` and `|C| = 14`, so
`T_AS->C = 1/(17*14) = 1/238 = 0.004201680672268907`, and `T_C->AS = 0`. Both
are the published values.

### 3.5 Derived combinations

```
T_inter^sum  = T_{A->B} + T_{B->A}
T_intra^sum  = T_A^mod  + T_B^mod
T_intra^avg  = (T_A^mod + T_B^mod) / 2
```

Per-module roll-ups in the module-level table:

```
T_module_out_sum(M) = sum over other modules O of T_{M->O}
T_module_in_sum(M)  = sum over other modules O of T_{O->M}
```

---

## 4. Impurity — `impurity()` — **56/56 exact**

Over the **core** subgraph (both endpoints team-labelled) with **self-loops
removed** — call that edge set `E`:

```
            #{ (i,j) in E : same team, A[i,j] < 0 }  +  #{ (i,j) in E : different teams, A[i,j] > 0 }
Impurity =  ------------------------------------------------------------------------------------------
                                              |E|
```

Edges are counted, not influence-weighted, so this is `l_max` independent.

## 5. Density — `density_core()`, `density_split()` — **56/56 exact each**

Same `E` and the same core node set `V_c`, `n = |V_c|`, team sizes `n_1, n_2`:

```
Density       = |E| / ( n (n-1) )

                #{ edges within a team }
intra_density = ------------------------
                n_1(n_1-1) + n_2(n_2-1)

                #{ edges between the teams }
inter_density = ----------------------------
                        2 n_1 n_2
```

**Departure — `N(N-1)`, not `N^2`.** The sheet writes `Density = total edges /
N^2`. The published numbers require `n(n-1)` on the core subgraph with
self-loops dropped. Worked check: `Comb_AS_C` has 25 core nodes and, after
dropping 1 self-loop, 64 core edges; `64/(25*24) = 0.10666666666666667`, the
published value, whereas `64/25^2 = 0.1024` and `75/31^2 = 0.0780` are not.
`impurity` uses the same 64-edge denominator: `5/64 = 0.078125`, also published.
Each denominator is the number of *ordered pairs that could carry an edge*, and
intra + inter recover the total: `38 + 26 = 64`.

---

## 6. State-based phenotypic score — `pheno_score_state_based()`

Per state `S`, per module, with the module's teams:

```
P_s^state(S) = mean_{i in T_1} S_i  -  mean_{i in T_2} S_i
```

Range `[-2, 2]`. No couplings enter; no normalisation by `N`.

**Departure — the +/-1 coding doubles the sheet's score.** The formula image
defines `E(T)` as *the fraction of nodes of T active in the state* and the
module score as `|E(A_1) - E(A_2)|`, which ranges over `[0, 1]` (signed,
`[-1, 1]`). `decode_state()` codes the state as `+/-1`, so a team mean is
`2 * E(T) - 1` and

```
mean(S over T_1) - mean(S over T_2)  =  2 * ( E(T_1) - E(T_2) )
```

exactly. Verified on every state of `Comb_RS_RS`: a fully committed state has
`E_1 - E_2 = 1.000` and `P_s = 2.000`, ratio 2.0 throughout.

**What this does and does not affect.** A global factor of 2 is invariant
under Pearson/Spearman correlation, so every correlation in
`Metrics/*.csv`, the `correlations_*` comparisons and the
`Figures/correlation_heatmaps/` panels is unchanged by it. It *does* matter
wherever an absolute cut or an absolute value is quoted: "fully committed"
is `|P_s| = 2` here but `1` on the sheet, and the hybrid-state thresholds
`0.5 / 1.0 / 1.5` used in `Figures/fig1_remaining.R` correspond to
`0.25 / 0.5 / 0.75` in the sheet's units. Divide by 2 to read any of these
in formula-sheet units.

**Note on the published column.** The published `_stateps` values are this
quantity divided by `N`: `dominant` reproduces as `published * N_network` for
170/180 module-network pairs. The extra `1/N` was dropped per instruction, so
the new numbers are `N` times the old ones.

The 10 non-reproducing pairs are in `Comb_C_E`, `Comb_C_M`, `Comb_E_RS`,
`Comb_M_PS`, `Comb_PS_RS`; all are exact sign flips, and in each the two most
populated states are near-tied mirror attractors (0.27662 vs 0.27638; 0.45273
vs 0.45234; 0.20016 vs 0.19884; 0.11036 vs 0.11011; 0.16535 vs 0.16385), so
which state counts as "dominant" is decided in the 4th decimal of a
Monte-Carlo basin estimate.

---

## 7. Influence-weighted phenotypic score — `teams_score()`, `pheno_score_state()` — **720/720 exact on `dominant`**

Per state `S`, per module, at a given `l_max`:

```
            1              1
F_T  =  ------- *  sum   ( --- *  sum   I[i,j] * S_j )
         |T|      i in T    N     j in V

P_s^infl(S)  =  | F_{T_1} - F_{T_2} |
```

Three points, each of which was wrong in my first attempt:

1. **The coupling is `I` at the given `l_max`, not `A`.** This is the whole
   meaning of "influence-weighted", and it is what makes the score `l_max`
   dependent. The published column varies with `l_max` in all 180
   (network, module) pairs — e.g. `Comb_AS_C`/AS runs 0.2488, 0.6485, 0.7745,
   0.8712 over `l_max` 1, 5, 7, 10 — so any `A`-based version is excluded
   immediately.
2. **The direction is outgoing.** `I[i,j]` with `i` in the team: the influence
   the team node *exerts*, weighted by the target's state. The incoming version
   (`sum_j I[j,i] S_j`, the field arriving at `i`) gives 0.2016 where the
   published value is 0.2488 at `l_max = 1`.
3. **The inner sum is a mean over all `N` network nodes**, i.e. the `1/N`. This
   was the missing factor of ~26. Note `j` ranges over the *whole network*, not
   just the module.

The absolute value in the last line is why the published column has no negative
entries (0 of 720) and why `dominant` reproduces 720/720 only once `|.|` is
applied. It also means this score is a *magnitude* and carries no direction
information, unlike §6.

### Range

`I[i,j] in [-1, 1]` and `S_j in {-1, +1}`, so `F_T in [-1, 1]` and
`P_s^infl in [0, 2]` — the same nominal span as §6, but reached only if every
influence were saturated. In practice the corpus spans 0.010 to 1.496 at the
module level and 0.101 to 1.619 at the network level. Worked example,
`Comb_RS_RS` at `l_max` 10: the committed state has `F_T1 = +0.9286`,
`F_T2 = -0.9095`, so `P_s^infl = 1.8381` where the state-based score is
exactly 2.00.

### The part to be careful about — this score is mostly structure

Unlike §6, which sees only the state, this score multiplies the state by the
influence matrix. For a **fully committed** state (`S = +1` on `T_1`, `-1` on
`T_2`) the state contributes nothing but signs, and

```
P_s^infl  =  (1/N) * |  |T_1|(M_11 - M_21)  -  |T_2|(M_12 - M_22)  |
```

where `M_kl` is the mean of `I` over block `(k,l)` — the *same four block
means* that `T_network` averages (§3.2). On a committed state the two
quantities are near-algebraically linked.

This is measurable. Evaluating `P_s^infl` on a **synthetic** committed state,
written down from the teams with no simulation at all:

| quantity | r against `T_network`, 56 networks, `l_max` 10 |
|---|---|
| `P_s^infl` on a synthetic committed state (pure structure) | **0.996** |
| `P_s^infl` on the actual simulated dominant state | 0.934 |
| `P_s^state` on the actual simulated dominant state | 0.044 |

The synthetic version, which contains no simulation output whatsoever,
correlates with `T_network` *better* than the simulated one does. And in
**30 of the 56 networks the dominant state simply is the fully committed
state**, so for those the two are the same number: the simulation adds
nothing. Restricting to the 26 networks where the dominant state is *not*
committed, the actual correlation falls to **0.564** while the synthetic one
holds at 0.975.

### What this score actually measures — the right reading

Calling the above "circularity" is too blunt. The matrix-weighted score is
**a state/topology alignment measure**, and it is exactly a bilinear form.
With `w` the *size-normalised* team vector (`w_i = +1/|T_1|` on team 1,
`-1/|T_2|` on team 2):

```
P_s^infl  =  | w^T I S | / N          (exact; verified to 2.8e-17)
```

and the frustration already stored in `frust0` is the same machinery with the
state in the left slot:

```
frustration  =  ( 1 - S^T A S / |E| ) / 2     (exact where both use the same node set)
```

Both project the vector `(matrix) * S` onto something; they differ only in
what:

| quantity | form | asks |
|---|---|---|
| frustration | `S^T A S` | does the field align with **the state**? |
| `P_s` (adjacency/influence) | `w^T A S` / `w^T I S` | does the field align with **the team partition**? |

They coincide exactly when `S = w`, i.e. on the fully committed state — which
is why 30 of 56 networks showed `P_s` equal to its synthetic-committed value,
and why that value tracks `T_network` so closely.

So the score is **not a corrupted landscape readout**; it is a team-structure
alignment measure, and being a joint function of topology and state is its
purpose rather than a defect. It is finer-grained than frustration: it is
team-resolved (`F_T1` and `F_T2` are separately meaningful), signed before the
`|.|`, and normalised per team so the larger team does not dominate. The
difference `F_T1 - F_T2` is what projects onto the team axis, on the
expectation that the two teams are driven oppositely.

**What this does and does not license.** Correlating `P_s^infl` with
`T_network` is a legitimate question — *do the attractors align with the team
partition, and more so when teams are stronger?* — and the answer is yes. What
it does **not** support is the looser reading that team strength predicts a
*landscape outcome*, because both sides are functions of the same topology:
controlling for the same score on a synthetic committed state leaves partial
r = +0.148. State it as an alignment result, not as a phenotype prediction.
The state-based score (§6) contains no matrix and remains the readout to use
when an outcome independent of topology is wanted; its r ~ 0.04 against
`T_network` is then the honest measurement, not a failure of the measure.

---

## 7b. Network-level phenotypic score — `network_ps_table()` — **80/80 exact on `dominant`**

Same shape as §7, but with the **network's own two teams** (from `<net>.teams`)
and evaluated entirely on the **core** network:

```
            1                1
F_T  =  --------- *  sum   ( ------ *  sum      I[i,j] * S_j )
         |T ∩ c|     i in     |c|     j in c
                     T ∩ c

P_s^net(S)  =  | F_{T_1} - F_{T_2} |
```

where `c` is the core node set — nodes surviving the iterative removal of
peripherals (§1) — and `|c| = n_core`.

**The core/full asymmetry is the thing to check here.** The module-level score
(§7) uses the *full* node set and `1/N`; this one uses the *core* set and
`1/n_core`. That is not a choice I would have predicted, but the full-node
version does not reproduce the published coefficients and the core version does.
For doubles at `l_max = 1` against the five published `dominant_network`
coefficients:

| variant | matches |
|---|---|
| **core matrix, `1/n_core`** | **5/5** |
| full matrix, `1/N` | 0/5 (0.9051 vs published 0.9502) |
| core node count but full matrix | 0/5 |
| no normalisation | 0/5 |
| mean of the module `|dF|` values | 0/5 |
| max of the module `|dF|` values | 0/5 |

Using `<net>.teams` is consistent with this: those files label only core nodes,
so the network partition is defined on `c` in the first place.

**How this was checked.** The per-network values were never saved — only the
`dominant_network` / `top30_network` / `top60_network` rows of
`correlations_pearson_equal_*.csv`, which are Pearson correlations of the
network-level P_s against the five network columns, computed within each module
count. So the check is indirect: recompute P_s per network, correlate, compare
coefficients. `dominant` reproduces **80/80** (4 module counts x 4 `l_max` x 5
columns; worst deviation 4.9e-05, the rounding in the published files).
`top30` and `top60` reproduce 0/80 each — the same basin-weight cause as §8.

`core_nodes()` was verified to reproduce `InfluenceMatrix()`'s pruning for all
56 networks, and `I[core, core]` to equal its reduced matrix to 0 difference,
so no extra call to that function (and none of its file side effects) is needed.

---

## 8. Aggregating over states — `aggregate_ps()`

States are sorted by basin size `f` descending; `cum_k` is the cumulative sum.
For a threshold `θ ∈ {0.30, 0.60}` let `K(θ)` be the smallest `k` with
`cum_k >= θ` (so the state that crosses the threshold **is included**).

```
dominant  =  P_s of the most populated state

                   sum_{k <= K(θ)} P_s,k * f_k
top-θ     =        ---------------------------      (basin-weighted mean)
                   sum_{k <= K(θ)} f_k
```

with `P_s` replaced by `|P_s|` for the influence-weighted variant.

Both variants use the basin-weighted **mean**; the weights are renormalised
over the selected states. They differ only in the absolute value:

| column | `absolute` | `normalise` |
|---|---|---|
| `P_s_state_*` | FALSE — signed, range [-2, 2] | TRUE |
| `P_s_infl_*` | TRUE — magnitude, `|dF|` | TRUE |

A weighted *sum* is excluded for both: a sum over states covering 30% of the
state space is roughly `0.3x` the typical `|P_s|` and would fall outside
`[min P_s, max P_s]` over the selected states, whereas the published value lies
inside that interval in every non-degenerate case (§ below).

### Why `top30` / `top60` do not reproduce exactly

For **both** variants the failures split perfectly on whether the basin weights
matter at all:

| variant | weights irrelevant (all selected states share one `P_s`) | weights matter |
|---|---|---|
| influence-weighted | 223/223 top30, 59/59 top60 reproduced | 0/497, 0/661 |
| state-based | 2/2 top30 reproduced | 0/178 top30, 0/180 top60 |

and in every non-reproducing case the published value still lies inside
`[min P_s, max P_s]` over the same selected states — 497/497 and 661/661
(influence-weighted), 178/178 and 180/180 (state-based). So the per-state
scores are right and only the *weights* differ: median relative error 1.0%
(top30) and 0.34% (top60). The `Avg0` values in the current `finFlagFreq.csv`
files are therefore not the ones the published columns were computed from —
ordinary Monte-Carlo scatter between Boolean runs. This is the same cause as
the 10 sign flips in §6.

---

## 9. States — `decode_state()`

`finFlagFreq.csv` gives one row per Boolean steady state: `states` is an
underscore-separated 0/1 string, `Avg0` the basin size (sums to 1 over the
file). The string is indexed by `<net>_nodes.txt`, in order; the length is
asserted to match. Coding into `S` is `1 -> +1`, `0 -> -1`.

---

## 10. Summary of checks

| quantity | check |
|---|---|
| `influence_full` vs `InfluenceMatrix` on the core block | max abs difference 0 |
| `density`, `impurity`, `intra_density`, `inter_density` | 56/56 networks each |
| `T_network` | 224/224 (56 networks x 4 `l_max`) |
| `T_module_intra` | 720/720 |
| `T_{A->B}` | 1800/1800 |
| influence-weighted `P_s` `dominant` (module level) | 720/720 |
| network-level `P_s` `dominant` (via published correlations) | 80/80 |
| state-based `P_s` `dominant` (vs published x `N`) | 170/180, the 10 being basin ties |
| `top30` / `top60` | reproduce exactly iff weight-independent (§8) |

Reproduce with `Rscript Metrics/validate_metrics.R`,
`validate_ps_inter.R`, `validate_infl_ps.R`, `validate_state_ps.R`,
`validate_network_ps.R`.
