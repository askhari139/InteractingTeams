## metrics_lib.R -----------------------------------------------------------
## Implements the metrics in Tkl_component_scatterplots/.../outline.md (the
## formula sheet): team strength, inter-module team strength, the inter/intra
## sums, teams score, phenotypic score, impurity and density.
##
## Conventions were pinned against the existing correlations_* outputs; see
## Metrics/README.md for which ones are confirmed and which are assumptions.
## -------------------------------------------------------------------------

suppressMessages(library(funcsKishore))

## Repo root. Override with options(interactingteams.repo = "/path") or by
## assigning REPO after sourcing.
REPO <- getOption("interactingteams.repo", getwd())

MODULE_DIR  <- "Single_EMT_MechSens"
NETWORK_DIR <- c("Double_EMT_MecnSens", "Triple_EMT_MechSens",
                 "Four_EMT_MechSens", "Five_EMT_MechSens")

## ---- IO ------------------------------------------------------------------

#' Read a .teams file: one comma-separated team per line.
read_teams <- function(path) {
  lapply(strsplit(readLines(path), ","), function(x) trimws(x[x != ""]))
}

#' Read a .topo file as a data frame with Type coded +1 / -1.
read_topo <- function(path) {
  df <- read.delim(path, sep = "", stringsAsFactors = FALSE)
  df$Type <- ifelse(df$Type == 2, -1L, 1L)
  df
}

#' Node order used by the Boolean state strings (the *_nodes.txt file).
read_nodes <- function(path) {
  n <- readLines(path)
  trimws(n[n != ""])
}

#' Locate a network file across the per-order directories.
find_net_file <- function(repo, net, ext) {
  for (d in NETWORK_DIR) {
    p <- file.path(repo, d, paste0(net, ".", ext))
    if (file.exists(p)) return(p)
  }
  stop("no ", ext, " file for ", net)
}

module_nodes <- function(repo, m) {
  read_nodes(file.path(repo, MODULE_DIR, paste0("Single_", m, "_nodes.txt")))
}

#' Module abbreviations encoded in a network name, e.g. Comb_AS_C_E.
modules_of <- function(net) {
  unique(strsplit(sub("^Comb_", "", net), "_")[[1]])
}

## ---- influence matrix ----------------------------------------------------

#' Influence matrix over *all* nodes, i.e. funcsKishore::InfluenceMatrix()
#' without the peripheral-pruning step.
#'
#' Reproduces InfluenceMatrix() exactly on the core sub-block (verified to
#' 0 difference), but keeps peripheral rows/cols, which the inter-module
#' metrics need: the single AS->C edge in Comb_AS_C lands on a peripheral.
influence_full <- function(topo_path, lmax) {
  im   <- TopoToIntMat(topo_path)
  mat  <- im[[1]]
  nds  <- im[[2]]
  amax <- abs(mat)
  res  <- 0
  for (l in seq_len(lmax)) {
    a <- ComputePowerMatrix(mat,  l)
    b <- ComputePowerMatrix(amax, l)
    r <- a / b
    r[is.nan(r)] <- a[is.nan(r)]
    res <- res + r
  }
  res <- res / lmax
  dimnames(res) <- list(nds, nds)
  res
}

## ---- team bookkeeping ----------------------------------------------------

#' Extend a partition so it covers `nodes`.
#'
#' The .teams files only label core nodes (peripherals are pruned before team
#' detection), but the intra-module team strength is computed over every module
#' node. An unlabelled node is assigned to the team it is most positively
#' coupled to, using the signed influence in both directions.
assign_unlabelled <- function(groups, nodes, I) {
  labelled <- unlist(groups)
  missing  <- setdiff(intersect(nodes, rownames(I)), labelled)
  for (p in missing) {
    score <- vapply(groups, function(g) {
      g <- intersect(g, rownames(I))
      if (!length(g)) return(-Inf)
      mean(I[p, g]) + mean(I[g, p])
    }, numeric(1))
    k <- which.max(score)
    groups[[k]] <- c(groups[[k]], p)
  }
  groups
}

#' Restrict a partition to a node set, dropping empty teams.
restrict_teams <- function(groups, nodes) {
  Filter(length, lapply(groups, function(g) intersect(g, nodes)))
}

## ---- team strength -------------------------------------------------------

#' T_kl block matrix.
#'
#' T_kl = | mean of I_ij over i in team k, j in team l |
#'
#' NOTE: the formula sheet writes T_kl as the mean of |I_ij|. The pipeline that
#' produced correlations_*/ instead takes the absolute value of the *block
#' mean*, matching funcsKishore::getGsVec(). We follow the implementation,
#' because that is what reproduces the published numbers; `abs_first = TRUE`
#' gives the literal formula-sheet reading.
Tkl <- function(I, groups, abs_first = FALSE) {
  k <- length(groups)
  out <- matrix(NA_real_, k, k, dimnames = list(names(groups), names(groups)))
  for (a in seq_len(k)) for (b in seq_len(k)) {
    blk <- I[groups[[a]], groups[[b]], drop = FALSE]
    out[a, b] <- if (abs_first) mean(abs(blk)) else abs(mean(blk))
  }
  out
}

#' T_network = (T00 + T01 + T10 + T11) / 4, over the network's own teams.
T_network <- function(I, net_teams, abs_first = FALSE) {
  g <- restrict_teams(net_teams, rownames(I))
  mean(Tkl(I, g, abs_first))
}

#' Intra-module team strength T_M^mod: the same 4-block average, evaluated on
#' the *network's* influence matrix but split by the *module's own* teams.
#'
#' Using the module's own partition (not the network's) is what makes T_M^mod
#' network-independent at lmax = 1, as the published outputs are: at lmax = 1
#' no path leaves the module, so the within-module block is the module's own
#' adjacency. At larger lmax the network matrix carries paths through the other
#' modules, so T_M^mod becomes network-dependent -- again as published.
T_module_intra <- function(I, mod_nodes, mod_teams, abs_first = FALSE) {
  g <- restrict_teams(mod_teams, mod_nodes)
  g <- restrict_teams(g, rownames(I))
  if (length(g) < 2) return(NA_real_)
  mean(Tkl(I, g, abs_first))
}

#' Node -> team lookup from allNodes.nodes (cached).
#'
#' allNodes.nodes is the authoritative module-level team definition: `Program`
#' is the module, `NewProgram` splits it into its two teams (e.g. Phase_Switch
#' -> Phase_Switch / Phase_acyclic). It covers every module node, peripherals
#' included, so no team-assignment heuristic is needed.
#'
#' Prefer this over the per-module `.teams` files: those disagree for PS
#' (6/9 vs the 12/3 here), and only the NewProgram split reproduces the
#' published PS numbers.
.newprogram_cache <- new.env(parent = emptyenv())
newprogram_map <- function(repo) {
  if (is.null(.newprogram_cache$map)) {
    nd <- read.csv(file.path(repo, "allNodes.nodes"), stringsAsFactors = FALSE)
    .newprogram_cache$map <- setNames(nd$NewProgram, nd$Node)
  }
  .newprogram_cache$map
}

#' Which NewProgram label is team 1, i.e. the positive pole of P_s.
#'
#' P_s = mean(S over team1) - mean(S over team2), so this table fixes the sign
#' of every phenotypic score. These orientations reproduce the published
#' `_stateps` signs: 30/30 for AS and 28/30 for C, M, PS and RS. EMT is
#' oriented mesenchymal-positive; taking EMT_Epi first flips the sign in 28 of
#' 30 networks.
TEAM1 <- c(Apoptotic_Switch   = "Apoptotic_Pro",
           CIP                = "CIP_E",
           EMT                = "EMT_Mes",
           Migration          = "Migration_Pro",
           Phase_Switch       = "Phase_Switch",
           Restriction_Switch = "Restriction_Cyclin")

#' The module's own 2-team partition, from allNodes.nodes.
#'
#' Returns a list of teams, team 1 first per TEAM1. Modules with a single
#' NewProgram (A, CC, GBM, GP, Gm, OL, Tbs) yield one team and therefore have
#' no two-team structure.
module_teams <- function(repo, m) {
  nodes <- module_nodes(repo, m)
  lab   <- unname(newprogram_map(repo)[nodes])
  g     <- split(nodes, lab)
  if (length(g) < 2) return(g)
  first <- intersect(TEAM1, names(g))
  if (length(first)) g <- g[c(first[1], setdiff(names(g), first[1]))]
  g
}

#' TRUE if the module has a genuine two-team split.
has_two_teams <- function(repo, m) length(module_teams(repo, m)) >= 2

#' Inter-module team strength T_{A->B}.
#'
#' T_{A->B} = | mean of I_ij over i in A, j in B |, on the network's full
#' influence matrix (peripherals kept -- in Comb_AS_C the only AS->C edge lands
#' on a peripheral node, so pruning would zero this out).
#'
#' As with Tkl, the absolute value is taken *after* the block mean, so opposing
#' cross-module edges cancel. `abs_first = TRUE` gives the literal
#' formula-sheet reading (mean of |I_ij|), which does not cancel.
T_inter <- function(I, nodes_a, nodes_b, abs_first = FALSE) {
  a <- intersect(nodes_a, rownames(I))
  b <- intersect(nodes_b, rownames(I))
  if (!length(a) || !length(b)) return(NA_real_)
  blk <- I[a, b, drop = FALSE]
  if (abs_first) mean(abs(blk)) else abs(mean(blk))
}

#' Team-aware inter-module strength.
#'
#' `T_inter()` above averages I over the WHOLE A x B block, so a coupling that
#' is perfectly aligned with the team structure -- A's team 1 driving B's team 1
#' and opposing B's team 2 -- has its two halves cancel and reports ~0. Over the
#' 450 ordered module pairs at l_max 10 the block-resolved value is larger in
#' 361 of them, median 4.3x; the extreme is M->C in Comb_C_M, 0.0714 against
#' 0.740. The two correlate only 0.61, so this is a reordering, not a rescaling.
#'
#' Both functions below resolve the four team sub-blocks first. Let
#' `m_kl = mean(I[A_k, B_l])`.
#'
#'   T_inter_blocks  = mean(|m_11|, |m_12|, |m_21|, |m_22|)
#'       magnitude of team-resolved coupling, the direct analogue of
#'       `T_network()` (§3.2), which averages the four |blocks| the same way.
#'
#'   T_inter_aligned = (m_11 - m_12 - m_21 + m_22) / 4
#'       SIGNED coherence with the team structure: positive when A's team 1
#'       drives B's team 1 and opposes B's team 2. This is exactly
#'       `w_A^T I w_B / 4` for the size-normalised team vectors
#'       `w_A = +1/|A_1| on A_1, -1/|A_2| on A_2` (and likewise `w_B`), so it
#'       is the inter-module member of the same bilinear family as the
#'       phenotypic scores (§6, §7).
#'
#' `|T_inter_aligned| == T_inter_blocks` exactly when all four sub-blocks carry
#' the aligned sign pattern, i.e. a perfectly coherent coupling.
#'
#' CAUTION on the sign. `w_A` depends on which team of A is called team 1, so
#' relabelling a module's teams flips the sign of `T_inter_aligned`. The
#' orientation is fixed by the `TEAM1` table above (EMT mesenchymal-positive,
#' etc.), so values are consistent across this project -- but the sign is only
#' interpretable relative to that convention. `|T_inter_aligned|` and
#' `T_inter_blocks` are convention-free; prefer them unless the orientation is
#' being asserted deliberately.
.inter_blocks <- function(I, teams_a, teams_b) {
  A1 <- intersect(teams_a[[1]], rownames(I)); A2 <- intersect(teams_a[[2]], rownames(I))
  B1 <- intersect(teams_b[[1]], colnames(I)); B2 <- intersect(teams_b[[2]], colnames(I))
  if (!length(A1) || !length(A2) || !length(B1) || !length(B2)) return(NULL)
  c(m11 = mean(I[A1, B1, drop = FALSE]), m12 = mean(I[A1, B2, drop = FALSE]),
    m21 = mean(I[A2, B1, drop = FALSE]), m22 = mean(I[A2, B2, drop = FALSE]))
}

T_inter_blocks <- function(I, teams_a, teams_b) {
  m <- .inter_blocks(I, teams_a, teams_b)
  if (is.null(m)) return(NA_real_)
  mean(abs(m))
}

T_inter_aligned <- function(I, teams_a, teams_b) {
  m <- .inter_blocks(I, teams_a, teams_b)
  if (is.null(m)) return(NA_real_)
  unname((m[["m11"]] - m[["m12"]] - m[["m21"]] + m[["m22"]]) / 4)
}

T_inter_sum  <- function(tab, a, b) tab[a, b] + tab[b, a]
T_intra_sum  <- function(vec, a, b) vec[[a]] + vec[[b]]
T_intra_mean <- function(vec, a, b) (vec[[a]] + vec[[b]]) / 2

## ---- impurity and density ------------------------------------------------

#' Edges of the core (team-labelled) subgraph, self-loops dropped.
core_edges <- function(topo, teams) {
  core <- unlist(teams)
  topo[topo$Source %in% core & topo$Target %in% core &
       topo$Source != topo$Target, , drop = FALSE]
}

#' Impurity = (negative edges within teams + positive edges between teams)
#'            / total edges
impurity <- function(topo, teams) {
  e   <- core_edges(topo, teams)
  lab <- setNames(rep(seq_along(teams), lengths(teams)), unlist(teams))
  same <- lab[e$Source] == lab[e$Target]
  sum((same & e$Type < 0) | (!same & e$Type > 0)) / nrow(e)
}

#' Density = total edges / N(N-1) over the core subgraph.
#'
#' The formula sheet writes N^2; the published numbers use the no-self-loop
#' denominator N(N-1), which is what we implement (56/56 networks reproduce).
density_core <- function(topo, teams) {
  e <- core_edges(topo, teams)
  n <- length(unlist(teams))
  nrow(e) / (n * (n - 1))
}

#' Within-team and between-team densities, same core subgraph.
density_split <- function(topo, teams) {
  e   <- core_edges(topo, teams)
  lab <- setNames(rep(seq_along(teams), lengths(teams)), unlist(teams))
  same <- lab[e$Source] == lab[e$Target]
  sz   <- lengths(teams)
  c(intra = sum(same)  / sum(sz * (sz - 1)),
    inter = sum(!same) / (2 * prod(sz)))
}

## ---- teams score and phenotypic score ------------------------------------

#' Teams score F_T, influence-weighted.
#'
#'   F_T = (1/|T|) * sum_{i in T} [ (1/N) * sum_{j in network} I_ij * S_j ]
#'
#' i.e. for each node i of the team, the mean influence it exerts on the
#' network weighted by each target's state, then averaged over the team.
#'
#' Three details, all needed to reproduce the published `influence_weighted`
#' column (verified at l_max 1 and 10 independently):
#'   * the coupling is the **influence matrix at the given l_max**, not the
#'     signed adjacency -- this is what makes the score l_max dependent;
#'   * the direction is **outgoing** (row i of I, i's influence on j), not the
#'     field arriving at i;
#'   * the inner sum is a **mean over all N network nodes**, not a bare sum.
teams_score <- function(I, state, team, n_network = ncol(I)) {
  team <- intersect(team, rownames(I))
  if (!length(team)) return(NA_real_)
  mean(as.vector(I[team, , drop = FALSE] %*% state[colnames(I)])) / n_network
}

#' Influence-weighted phenotypic score of one state: dF = F_T1 - F_T2.
pheno_score_state <- function(I, state, teams, n_network = ncol(I)) {
  teams_score(I, state, teams[[1]], n_network) -
  teams_score(I, state, teams[[2]], n_network)
}

#' Decode a finFlagFreq `states` string ("1_0_1_...") into a named +/-1 vector.
decode_state <- function(s, nodes) {
  bits <- as.integer(strsplit(s, "_")[[1]])
  stopifnot(length(bits) == length(nodes))
  setNames(ifelse(bits == 1, 1, -1), nodes)
}

#' Per-state influence-weighted scores with basin sizes.
#'
#' `I` must be the network influence matrix at the l_max of interest, so this
#' table has to be rebuilt for each l_max (unlike the state-based scores).
#' Returns a data frame: state, freq (Avg0, basin size), ps.
state_pheno_scores <- function(freq_path, I, nodes, teams) {
  fr <- read.csv(freq_path, stringsAsFactors = FALSE)
  Is <- I[nodes, nodes, drop = FALSE]
  ps <- vapply(fr$states, function(s)
    pheno_score_state(Is, decode_state(s, nodes), teams, length(nodes)),
    numeric(1))
  data.frame(state = fr$states, freq = fr$Avg0, ps = as.numeric(ps),
             stringsAsFactors = FALSE)
}

#' Aggregate per-state scores the three ways the pipeline uses.
#'
#'  dominant : P_s of the single most populated state (PStop)
#'  top30    : aggregate over the states making up the top 30% of the state
#'             space, ranked by basin size (PSw)
#'  top60    : same, 60%
#'
#' @param absolute  take |P_s| per state before aggregating. The published
#'   influence-weighted column is |dF| (0 of 720 values are negative, and
#'   `dominant` reproduces 720/720 once the absolute value is taken).
#' @param normalise TRUE renormalises the basin weights over the selected
#'   states (a weighted *mean*); FALSE leaves them unnormalised (a weighted
#'   *sum*, which is the stated definition for the state-based score).
#'
#' The state that crosses the threshold is included.
aggregate_ps <- function(df, absolute = FALSE, normalise = FALSE) {
  if (absolute) df$ps <- abs(df$ps)
  d   <- df[order(-df$freq), ]
  cum <- cumsum(d$freq)
  pick <- function(thresh) {
    k <- which(cum >= thresh)[1]
    if (is.na(k)) k <- nrow(d)
    sub <- d[seq_len(k), ]
    if (normalise) sum(sub$ps * sub$freq) / sum(sub$freq) else sum(sub$ps * sub$freq)
  }
  c(dominant = d$ps[1], top30 = pick(0.30), top60 = pick(0.60))
}

## ---- state-based phenotypic score ---------------------------------------

#' State-based phenotypic score of one state for one module.
#'
#' P_s = mean(S over team 1) - mean(S over team 2)
#'
#' The `_stateps` / "StateBased_newnorm" variant: it uses only the state vector,
#' no couplings. Already a mean difference, so no further normalisation by
#' network size. Range is [-2, 2] in +/-1 coding.
pheno_score_state_based <- function(state, teams) {
  t1 <- intersect(teams[[1]], names(state))
  t2 <- intersect(teams[[2]], names(state))
  if (!length(t1) || !length(t2)) return(NA_real_)
  mean(state[t1]) - mean(state[t2])
}

#' Per-state state-based scores with basin sizes, from a finFlagFreq file.
state_ps_table <- function(freq_path, nodes, teams) {
  fr <- read.csv(freq_path, stringsAsFactors = FALSE)
  ps <- vapply(fr$states, function(s)
    pheno_score_state_based(decode_state(s, nodes), teams), numeric(1))
  data.frame(state = fr$states, freq = fr$Avg0, ps = as.numeric(ps),
             stringsAsFactors = FALSE)
}

## ---- network-level phenotypic score -------------------------------------

#' The core node set: iteratively drop nodes that are not both a source and a
#' target. Replicates the pruning inside funcsKishore::InfluenceMatrix(), so
#' `I[core_nodes(topo), core_nodes(topo)]` equals that function's return value
#' (verified to 0 difference) without its file/`setwd` side effects.
core_nodes <- function(topo) {
  df <- topo
  repeat {
    nodes <- unique(c(df$Source, df$Target))
    peri  <- nodes[!(nodes %in% df$Source & nodes %in% df$Target)]
    if (!length(peri)) break
    df <- df[!(df$Source %in% peri | df$Target %in% peri), , drop = FALSE]
  }
  sort(unique(c(df$Source, df$Target)))
}

#' Network-level phenotypic score, per state.
#'
#' Same shape as the module-level influence-weighted score, but evaluated
#' entirely on the **core** (peripheral-pruned) network and using the
#' network's own two teams:
#'
#'   F_T = (1/|T|) sum_{i in T} [ (1/n_core) sum_{j in core} I_ij S_j ]
#'   P_s = | F_T1 - F_T2 |
#'
#' Note the asymmetry with the module-level score (§7 of FORMULAE.md), which
#' uses the *full* node set and `1/N`. Using the full matrix here does not
#' reproduce the published `dominant_network` correlations; the core version
#' reproduces all 80 of them exactly.
network_ps_table <- function(freq_path, I, nodes, net_teams, core) {
  fr <- read.csv(freq_path, stringsAsFactors = FALSE)
  core <- intersect(core, rownames(I))
  Ic   <- I[core, core, drop = FALSE]
  g    <- restrict_teams(net_teams, core)
  if (length(g) < 2) return(NULL)
  ps <- vapply(fr$states, function(s) {
    Sv <- decode_state(s, nodes)[core]
    f  <- function(T) mean(as.vector(Ic[T, , drop = FALSE] %*% Sv)) / length(core)
    abs(f(g[[1]]) - f(g[[2]]))
  }, numeric(1))
  data.frame(state = fr$states, freq = fr$Avg0, ps = as.numeric(ps),
             stringsAsFactors = FALSE)
}

## ---- PCA metrics ---------------------------------------------------------

#' PCA spectrum metrics for one network, from <Tag>_PCA_variance.csv.
#'
#' Definitions verified to reproduce the published correlations_*.csv
#' coefficients for PC1_var, NumPC_90, PCA_entropy and PCA_entropy_norm.
#'   PC1_var          variance explained by PC1
#'   NumPC_90         PCs needed to reach 90% cumulative variance
#'   PCA_entropy      Shannon entropy of the normalised spectrum
#'   PCA_entropy_norm PCA_entropy / log(n_PCs)
.pca_cache <- new.env(parent = emptyenv())
pca_metrics <- function(repo, net) {
  if (is.null(.pca_cache$tab)) {
    fs <- Sys.glob(file.path(repo, "*_EMT_Mec?Sens", "*_PCA_variance.csv"))
    tab <- list()
    for (f in fs) {
      d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      vc <- grep("^PC\\d+_variance_explained$", names(d), value = TRUE)
      vc <- vc[order(as.integer(sub("^PC(\\d+)_.*", "\\1", vc)))]
      for (i in seq_len(nrow(d))) {
        v <- suppressWarnings(as.numeric(unlist(d[i, vc])))
        v <- v[!is.na(v)]
        if (!length(v) || sum(v) <= 0) next
        p   <- v / sum(v)
        cum <- cumsum(p)
        H   <- -sum(p[p > 0] * log(p[p > 0]))
        tab[[d$topology[i]]] <- list(
          PC1_var = v[1], NumPC_90 = which(cum >= 0.90)[1],
          PCA_entropy = H,
          PCA_entropy_norm = if (length(v) > 1) H / log(length(v)) else NA_real_,
          n_PCs = length(v), n_samples = d$n_samples[i])
      }
    }
    .pca_cache$tab <- tab
  }
  .pca_cache$tab[[net]]
}

#' Network-level state-based phenotypic score, per state.
#'
#'   P_s = mean(S over team 1) - mean(S over team 2)
#'
#' The network analogue of pheno_score_state_based(): the network's own two
#' teams instead of a module's. Signed (range [-2, 2]) and l_max free, since no
#' couplings enter. `<net>.teams` labels only core nodes, so only core nodes
#' contribute -- the same restriction network_ps_table() makes explicit.
#'
#' Note the deliberate asymmetry with network_ps_table(), which returns |dF|:
#' the absolute value there is what reproduces the published `*_network`
#' coefficients, whereas the state-based score is kept signed to match its
#' module-level counterpart.
network_state_ps_table <- function(freq_path, nodes, net_teams) {
  fr <- read.csv(freq_path, stringsAsFactors = FALSE)
  g  <- restrict_teams(net_teams, nodes)
  if (length(g) < 2) return(NULL)
  ps <- vapply(fr$states, function(s)
    pheno_score_state_based(decode_state(s, nodes), g), numeric(1))
  data.frame(state = fr$states, freq = fr$Avg0, ps = as.numeric(ps),
             stringsAsFactors = FALSE)
}
