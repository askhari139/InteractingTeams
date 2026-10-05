# Clustering of networks by phenotypic score (Q3, Q4)

Answers Q3 and Q4 of
`Tkl_component_scatterplots/StateBased_newnorm/lmax10/analysis.md`.

Scripts: `clustering_analysis.R`, `clustering_plots.R`.
Outputs: `clustering_inconsistency.csv`, `clustering_null.csv`,
`clustering_null_B.csv`, `clustering_vs_properties.csv`,
`figures/clustering_*.png`.

## Method

The document specifies a method for Q1/Q2 (line 10) but not for Q3, and
"most inconsistent clustering" admits more than one reading. Three are computed
and kept separate:

| reading | question it answers | statistic |
|---|---|---|
| **A** stability | does removing this module change the clustering? | 1 − ARI between the full-corpus labels restricted to the 26 retained networks and the labels recomputed on those 26 |
| **B** coherence | do this module's networks group together? | normalised Shannon entropy of the cluster distribution of the 30 networks containing it (0 = one cluster, 1 = uniform) |
| **C** explanatory power | does membership explain the clustering? | ARI between the clustering and the binary "contains X" split |

Features: the eight P_s columns per network at a given l_max (module-level and
network-level × influence-weighted and state-based × top30/top60), z-scored.
Ward linkage on Euclidean distance, cut at k = 2, 3, 4. Both A and B carry a
random-subset null (400 draws), because neither statistic is interpretable
without one: A removes 30 of 56 networks, and B can look coherent purely
because the cluster sizes are uneven.

Cluster quality is moderate (mean silhouette 0.26–0.55) and the cluster sizes
are lopsided — e.g. 18/38 at l_max 10, k = 2 — which is why the nulls matter.

## Q3 — which module has the most inconsistent clustering?

**The three readings disagree, and two of them find nothing.**

### Reading A — stability: nothing to report

| module | mean 1 − ARI |
|---|---|
| PS | 0.638 |
| E | 0.486 |
| C | 0.422 |
| M | 0.414 |
| RS | 0.363 |
| AS | 0.362 |

The random-removal null has median 0.431 and a 5th–95th range of
**[0.000, 0.830]** — almost the whole scale. Only **3 of 72** module removals
beat the null 95th percentile, where chance gives ~4. The clustering is simply
unstable to dropping 30 of 56 networks, whichever 30. This reading cannot
separate the modules.

### Reading B — coherence: no module is scattered more than chance

| module | mean normalised entropy | largest cluster share |
|---|---|---|
| M | 0.723 | 0.669 |
| C | 0.681 | 0.697 |
| RS | 0.678 | 0.669 |
| E | 0.614 | 0.744 |
| AS | 0.476 | 0.794 |
| PS | **0.259** | **0.889** |

Null: median 0.717, 5th–95th [0.418, 0.935]. So M, C and RS — the apparent
"most inconsistent" — sit exactly at the null median. **0 of 72 module spreads
exceed the null 95th percentile**, while **25 of 72 fall below its 5th
percentile** (chance ~4 each). The signal runs the opposite way to the
question: no module clusters worse than chance, and several cluster
*better* than chance. PS is the extreme case — entropy 0.259 against a null of
0.717, and at l_max 10, k = 3 all 30 of its networks land in one cluster.

### Reading C — explanatory power: PS, and only PS

ARI between the clustering and "contains module X":

| module | l_max 1 | 5 | 7 | 10 |
|---|---|---|---|---|
| PS | 0.075 | 0.164 | 0.171 | **0.452** |
| AS | 0.075 | 0.109 | 0.062 | 0.009 |
| E | 0.040 | 0.005 | 0.007 | 0.004 |
| C | 0.057 | −0.003 | −0.015 | −0.012 |
| M | 0.015 | −0.016 | −0.013 | −0.016 |
| RS | −0.018 | −0.009 | −0.010 | −0.011 |

**PS is the only module whose membership tracks the clustering, and it
strengthens sharply with l_max — ARI 0.45 at l_max 10 against ~0 for all five
others.** At l_max 10 the primary split of networks in phenotypic-score space
is substantially "contains PS or not".

I checked the obvious over-reading: at l_max 1, k = 3 one cluster happens to
hold exactly 30 networks, but it is *not* the PS set (21 PS + 9 non-PS), so the
l_max 1 correspondence is weak, as the ARI of 0.075 says.

### The answer

Read as "whose networks fail to group" the nominal answer is **M**, but it is
indistinguishable from chance and should not be reported. The defensible finding
is the reverse: **PS is the module that most determines the clustering** — most
coherent by reading B, and the only one with explanatory power by reading C.
This is consistent with PS being the structural outlier found earlier: by far
the weakest team (T_s 0.086 against 0.50–0.92) and the only module with non-zero
internal impurity (0.143).

It also, finally, justifies the `_noPS` folders. PS barely affects the
*correlations* (Q1: joint-lowest impact), but it dominates the *clustering*.
Separating PS-containing networks is well motivated for landscape structure,
not for correlation strength.

## Q4 — which property correlates with the inconsistency?

Mean Spearman ρ over l_max and k, n = 6 modules per correlation:

| property | A: stability | B: coherence | C: explanatory power |
|---|---|---|---|
| impurity across modules | −0.27 | **+0.61** | **−0.51** |
| intra-team density | +0.04 | +0.42 | −0.48 |
| density within | −0.03 | +0.34 | −0.45 |
| density across modules | −0.09 | +0.38 | −0.44 |
| team strength | −0.13 | +0.20 | −0.31 |
| inter-team density | +0.14 | −0.22 | +0.09 |

**Cross-module impurity is the best candidate**, and it is the only property
with a coherent story across the two readings that carry signal: modules whose
inter-module edges are team-misaligned have their networks *scattered* across
clusters (B, ρ = +0.61) and *fail* to explain the clustering (C, ρ = −0.51).
The mechanism is plausible — misaligned coupling scrambles a module's
phenotypic score inconsistently from network to network.

The supporting numbers line up: C (impurity across 0.957) and M (0.842) are the
two most scattered modules; PS (0.357) is the most coherent. RS is the
exception — low impurity (0.360) but high scatter (0.678) — which is why ρ is
0.61 rather than near 1.

**But n = 6.** ρ = 0.61 is not significant at six points, `density_within` is a
deterministic weighted average of the two team densities so those three rows are
not independent, and reading A — the one that follows the document's own
leave-one-out method — shows nothing for any property. Treat cross-module
impurity as **the hypothesis worth testing**, not a result. Testing it needs
more modules (see `loo_results.md`, "What would make this answerable") or an
independent corpus.

## Figures

| file | content |
|---|---|
| `clustering_stability.png` | Reading A: 1 − ARI per module against the random-removal null violin, faceted by k. |
| `clustering_coherence.png` | Reading B: normalised entropy of cluster membership per module, faceted by k. |
| `clustering_membership_ARI.png` | Reading C: ARI of clustering vs contains-module, by l_max, one line per module. |
| `clustering_composition.png` | How each module's 30 networks distribute over 3 clusters at l_max 10, with counts. |
| `clustering_vs_properties.png` | All three readings against each module property, coloured by module, ρ annotated. |

## Caveat on the feature set

The eight P_s columns are far from independent: the influence-weighted pair
correlates ~0.93 with network T_s while the state-based pair correlates ~0, so
the clustering is driven mostly by the influence-weighted columns. A version
restricted to the two influence-weighted network-level columns, or one that
weights the families equally, could give a different clustering. That choice was
not specified in the document and is worth confirming before building on these
results.
