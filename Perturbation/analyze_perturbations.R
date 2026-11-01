## analyze_perturbations.R ------------------------------------------------
## Step 4 of the perturbation study on Comb_AS_C_E_M_RS: from the Boolean runs
## (<net>_finFlagFreq.csv + <net>_nodes.txt) and the topo files, compute for the
## original and all 30 perturbed networks
##
##   structure   network T_s, intra-module T_s (per module, mean, sum),
##               inter-module T_s (|block mean| = the pipeline's T_A->B, and the
##               team-resolved version), their sums/means, density, impurity
##   landscape   phenotypic score P_s: module-level (influence-weighted and
##               state-based, top 30% and dominant) and network-level
##   PCA         PC1 variance, PCs to 90%, entropy -- ONLY if RACIPE
##               *_PCA_variance.csv files are found (see PCA_FILES below)
##
## and then draw scatter plots of the landscape against the structure and the
## correlation matrix, plus the three kinds against the original.
##
## Definitions are those of metrics_lib.R / FORMULAE.md, re-implemented here so
## the script needs only ggplot2 and cowplot (no funcsKishore):
##   influence matrix   I = (1/l_max) sum_l A^l / |A|^l  (0/0 -> 0), [source, target]
##   T_kl               |mean of I over team k x team l|   (abs AFTER the mean)
##   T_network          mean of the four T_kl over the network's own teams
##                      (the ORIGINAL network's .teams partition, held fixed)
##   T_module_intra     same, on the module's own teams (allNodes.nodes NewProgram)
##   T_A->B             |mean of I over module A x module B|, full (unpruned) matrix
##   P_s infl (module)  |F_T1 - F_T2| per state, F_T = (1/|T|) sum_{i in T}
##                      (1/N) sum_j I_ij S_j, module teams, N = all nodes;
##                      aggregated as a basin-weighted MEAN over the states that
##                      make up the top 30% of the state space
##   P_s state (module) mean(S over team 1) - mean(S over team 2), signed, summed
##                      with basin weights over the top 30% (unnormalised sum)
##   network-level      same two scores with the network's own teams, on the
##                      core (peripheral-pruned) network for the influence version
##
## CAUTION (FORMULAE.md s7): the influence-weighted P_s multiplies the state by
## the same influence matrix that defines T_s, so its correlation with
## structure is partly built in. The state-based score has no matrix in it and
## is the independent check. Both are drawn.
## CAUTION: one Boolean run per network. Basin sizes drift between runs of the
## same network (README, known gap 3), so differences between networks of the
## same kind are partly simulation noise.
##
## Validation: if Metrics/metrics_network_level.csv etc. exist and contain
## Comb_AS_C_E_M_RS, the original network's values are printed next to them.
##
## Run from the repo root:
##   Rscript Perturbation/analyze_perturbations.R
## or with explicit paths:
##   Rscript Perturbation/analyze_perturbations.R <topo_dir> <design_dir> <allNodes.nodes> <orig.teams>
##
## Outputs in Perturbation/analysis_<modules>/ :
##   perturbation_network_level.csv   perturbation_module_level.csv
##   perturbation_pair_level.csv      validation_original.csv (if the Metrics tables exist)
##   scatter_landscape_vs_structure.png  correlation_matrix.png  kinds_vs_original.png
## -------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(cowplot)})

args <- commandArgs(trailingOnly = TRUE)
topodir   <- if (length(args) >= 1) args[1] else file.path("Perturbation", "topo_AS_C_E_M_RS")
designdir <- if (length(args) >= 2) args[2] else file.path("Perturbation", "design_AS_C_E_M_RS")
nodes_f   <- if (length(args) >= 3) args[3] else "allNodes.nodes"
teams_f   <- if (length(args) >= 4) args[4] else file.path("Five_EMT_MechSens", "Comb_AS_C_E_M_RS.teams")
LMAX <- 10
PCA_FILES <- Sys.glob(file.path(topodir, "*_PCA_variance.csv"))   # RACIPE output, if you have it
mods <- c("AS", "C", "E", "M", "RS")
tag  <- paste(mods, collapse = "_")
outdir <- file.path("Perturbation", paste0("analysis_", tag)); dir.create(outdir, FALSE, TRUE)

PROG2MOD <- c(Apoptotic_Switch = "AS", CIP = "C", EMT = "E", Migration = "M",
              Phase_Switch = "PS", Restriction_Switch = "RS")
TEAM1 <- c(Apoptotic_Switch = "Apoptotic_Pro", CIP = "CIP_E", EMT = "EMT_Mes",
           Migration = "Migration_Pro", Phase_Switch = "Phase_Switch",
           Restriction_Switch = "Restriction_Cyclin")

## ---- helpers (same definitions as metrics_lib.R) -----------------------------------------
restrict_teams <- function(groups, nodes) Filter(length, lapply(groups, function(g) intersect(g, nodes)))
Tkl_mean <- function(I, groups) {
  k <- length(groups); if (k < 2) return(NA_real_)
  mean(sapply(seq_len(k), function(a) sapply(seq_len(k), function(b) abs(mean(I[groups[[a]], groups[[b]], drop = FALSE])))))
}
influence <- function(A, lmax) {
  Ab <- abs(A); P <- diag(nrow(A)); Pb <- diag(nrow(A)); res <- 0
  for (l in seq_len(lmax)) {
    P <- P %*% A; Pb <- Pb %*% Ab
    r <- P / Pb; r[is.nan(r)] <- P[is.nan(r)]
    res <- res + r
  }
  res <- res / lmax; dimnames(res) <- dimnames(A); res
}
core_nodes <- function(topo) {       # drop nodes that are not both a source and a target, iteratively
  df <- topo
  repeat {
    nodes <- unique(c(df$Source, df$Target))
    peri <- nodes[!(nodes %in% df$Source & nodes %in% df$Target)]
    if (!length(peri)) break
    df <- df[!(df$Source %in% peri | df$Target %in% peri), , drop = FALSE]
  }
  sort(unique(c(df$Source, df$Target)))
}
decode_matrix <- function(states, nodes) {        # states "1_0_1..." -> n_states x n matrix of +/-1
  bits <- do.call(rbind, strsplit(states, "_", fixed = TRUE))
  if (ncol(bits) != length(nodes)) stop("state strings have ", ncol(bits), " entries but there are ", length(nodes), " nodes")
  S <- matrix(ifelse(as.integer(bits) == 1L, 1, -1), nrow = nrow(bits)); colnames(S) <- nodes; S
}
agg_ps <- function(ps, freq, absolute, normalise) {   # = aggregate_ps() in metrics_lib.R
  if (absolute) ps <- abs(ps)
  o <- order(-freq); ps <- ps[o]; f <- freq[o]; cum <- cumsum(f)
  k <- which(cum >= 0.30)[1]; if (is.na(k)) k <- length(f)
  top30 <- if (normalise) sum(ps[1:k] * f[1:k]) / sum(f[1:k]) else sum(ps[1:k] * f[1:k])
  c(dominant = ps[1], top30 = top30)
}
## F_T for every state at once: (1/|T|) sum_{i in T} (1/N) sum_j I_ij S_j
team_score <- function(I, S, team, N) {
  team <- intersect(team, rownames(I)); if (!length(team)) return(rep(NA_real_, nrow(S)))
  as.vector(colMeans(I[team, colnames(S), drop = FALSE] %*% t(S))) / N
}

## ---- inputs ---------------------------------------------------------------------------------------
man <- read.csv(file.path(designdir, "manifest.csv"), stringsAsFactors = FALSE)
nd  <- read.csv(nodes_f, stringsAsFactors = FALSE)
net_teams <- lapply(strsplit(readLines(teams_f), ","), function(x) trimws(x[x != ""]))
cat(sprintf("%d networks in the manifest; topo dir %s\n", nrow(man), topodir))

nodeinfo <- function(nodes) {
  d <- data.frame(node = nodes, stringsAsFactors = FALSE)
  d$program <- nd$Program[match(d$node, nd$Node)]; d$team <- nd$NewProgram[match(d$node, nd$Node)]
  d$module <- unname(PROG2MOD[d$program])
  if (anyNA(d$module) || any(!d$module %in% mods)) stop("nodes outside ", tag, ": ",
        paste(d$node[is.na(d$module) | !d$module %in% mods], collapse = ", "))
  t1 <- TEAM1[d$program]; d$trank <- ifelse(!is.na(t1) & d$team == t1, 0L, 1L)
  d[order(match(d$module, mods), d$trank, d$team, d$node), ]
}
mod_teams_of <- function(info) lapply(setNames(mods, mods), function(m) {
  d <- info[info$module == m, ]
  split(d$node, factor(d$team, levels = unique(d$team[order(d$trank, d$team)])))
})

## ---- one network -----------------------------------------------------------------------------------
analyse <- function(net) {
  topo <- read.delim(file.path(topodir, paste0(net, ".topo")), sep = "", stringsAsFactors = FALSE)
  topo$sign <- ifelse(topo$Type == 2, -1, 1)
  state_nodes <- {x <- readLines(file.path(topodir, paste0(net, "_nodes.txt"))); trimws(x[x != ""])}
  fr <- read.csv(file.path(topodir, paste0(net, "_finFlagFreq.csv")), stringsAsFactors = FALSE)
  if (!all(c("states", "Avg0") %in% names(fr))) stop(net, " finFlagFreq needs columns states, Avg0; found: ", paste(names(fr), collapse = ", "))
  info <- nodeinfo(unique(c(topo$Source, topo$Target))); nodes <- info$node; n <- length(nodes)
  if (!setequal(state_nodes, nodes)) stop(net, ": nodes.txt and topo disagree on the node set")
  A <- matrix(0, n, n, dimnames = list(nodes, nodes))
  for (k in seq_len(nrow(topo))) A[topo$Source[k], topo$Target[k]] <- topo$sign[k]
  I <- influence(A, LMAX)
  mt <- mod_teams_of(info); idx <- lapply(setNames(mods, mods), function(m) info$node[info$module == m])
  S <- decode_matrix(fr$states, state_nodes); freq <- fr$Avg0
  if (abs(sum(freq) - 1) > 0.05) message(sprintf("  note: %s basin sizes sum to %.3f, not 1", net, sum(freq)))

  ## structure
  Tnet <- Tkl_mean(I, restrict_teams(net_teams, rownames(I)))
  intra <- sapply(mods, function(m) Tkl_mean(I, restrict_teams(mt[[m]], rownames(I))))
  pairs <- do.call(rbind, lapply(mods, function(a) lapply(mods, function(b) if (a != b) {
    ma <- mt[[a]]; mb <- mt[[b]]
    m4 <- c(mean(I[ma[[1]], mb[[1]], drop = FALSE]), mean(I[ma[[1]], mb[[2]], drop = FALSE]),
            mean(I[ma[[2]], mb[[1]], drop = FALSE]), mean(I[ma[[2]], mb[[2]], drop = FALSE]))
    data.frame(network = net, module_from = a, module_to = b,
               T_from_to = abs(mean(I[idx[[a]], idx[[b]], drop = FALSE])),
               T_blocks = mean(abs(m4)), T_aligned = (m4[1] - m4[2] - m4[3] + m4[4]) / 4)
  }) |> (\(l) Filter(Negate(is.null), l))()) |> (\(l) do.call(rbind, l))())
  out_sum <- tapply(pairs$T_from_to, pairs$module_from, sum)[mods]
  in_sum  <- tapply(pairs$T_from_to, pairs$module_to, sum)[mods]
  core_e <- topo[topo$Source %in% unlist(net_teams) & topo$Target %in% unlist(net_teams) & topo$Source != topo$Target, ]
  lab <- setNames(rep(seq_along(net_teams), lengths(net_teams)), unlist(net_teams))
  same <- lab[core_e$Source] == lab[core_e$Target]
  nc <- length(unlist(net_teams))
  dens <- nrow(core_e) / (nc * (nc - 1))
  imp  <- sum((same & core_e$sign < 0) | (!same & core_e$sign > 0)) / nrow(core_e)

  ## module-level P_s (module teams, all N nodes)
  modrows <- do.call(rbind, lapply(mods, function(m) {
    g <- mt[[m]]; N <- n
    dF <- team_score(I, S, g[[1]], N) - team_score(I, S, g[[2]], N)
    St <- rowMeans(S[, intersect(g[[1]], colnames(S)), drop = FALSE]) - rowMeans(S[, intersect(g[[2]], colnames(S)), drop = FALSE])
    ai <- agg_ps(dF, freq, absolute = TRUE,  normalise = TRUE)
    as <- agg_ps(St, freq, absolute = FALSE, normalise = FALSE)
    data.frame(network = net, module = m, T_module_intra = intra[[m]],
               out_sum = out_sum[[m]], in_sum = in_sum[[m]],
               ps_infl_top30 = ai[["top30"]], ps_infl_dominant = ai[["dominant"]],
               ps_state_top30 = as[["top30"]], ps_state_dominant = as[["dominant"]])
  }))
  ## network-level P_s (network's own teams)
  core <- intersect(core_nodes(topo), rownames(I)); g <- restrict_teams(net_teams, core)
  Sc <- S[, core, drop = FALSE]; Ic <- I[core, core, drop = FALSE]
  f <- function(T) as.vector(colMeans(Ic[T, , drop = FALSE] %*% t(Sc))) / length(core)
  ni <- agg_ps(abs(f(g[[1]]) - f(g[[2]])), freq, absolute = FALSE, normalise = TRUE)
  gs <- restrict_teams(net_teams, nodes)
  Stn <- rowMeans(S[, gs[[1]], drop = FALSE]) - rowMeans(S[, gs[[2]], drop = FALSE])
  ns <- agg_ps(Stn, freq, absolute = FALSE, normalise = FALSE)

  netrow <- data.frame(network = net, T_network = Tnet,
    T_inter_sum = sum(pairs$T_from_to), T_inter_mean = mean(pairs$T_from_to),
    T_inter_blocks_mean = mean(pairs$T_blocks), T_inter_blocks_sum = sum(pairs$T_blocks),
    T_intra_mean = mean(intra), T_intra_sum = sum(intra),
    density = dens, impurity = imp,
    ps_net_infl_top30 = ni[["top30"]], ps_net_infl_dominant = ni[["dominant"]],
    ps_net_state_top30 = ns[["top30"]], ps_net_state_dominant = ns[["dominant"]],
    ps_mod_infl_top30 = mean(modrows$ps_infl_top30), ps_mod_state_top30 = mean(modrows$ps_state_top30),
    ps_mod_state_abs_top30 = mean(abs(modrows$ps_state_top30)),
    n_states = nrow(fr))
  list(net = netrow, mod = modrows, pair = pairs)
}

res <- lapply(man$network, function(nm) { cat(" ", nm, "\n"); analyse(nm) })
NET  <- merge(man[, c("network", "kind", "rep", "n_removed", "n_added")], do.call(rbind, lapply(res, `[[`, "net")), by = "network")
NET  <- NET[order(match(NET$network, man$network)), ]
MODL <- do.call(rbind, lapply(res, `[[`, "mod"))
PAIR <- do.call(rbind, lapply(res, `[[`, "pair"))

## optional PCA metrics from RACIPE
if (length(PCA_FILES)) {
  tab <- list()
  for (f in PCA_FILES) {
    d <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    vc <- grep("^PC\\d+_variance_explained$", names(d), value = TRUE)
    vc <- vc[order(as.integer(sub("^PC(\\d+)_.*", "\\1", vc)))]
    for (i in seq_len(nrow(d))) {
      v <- suppressWarnings(as.numeric(unlist(d[i, vc]))); v <- v[!is.na(v)]
      if (!length(v) || sum(v) <= 0) next
      p <- v / sum(v); H <- -sum(p[p > 0] * log(p[p > 0]))
      tab[[d$topology[i]]] <- c(PC1_var = v[1], NumPC_90 = which(cumsum(p) >= 0.9)[1],
                                PCA_entropy_norm = if (length(v) > 1) H / log(length(v)) else NA_real_)
    }
  }
  for (cn in c("PC1_var", "NumPC_90", "PCA_entropy_norm"))
    NET[[cn]] <- vapply(NET$network, function(x) if (is.null(tab[[x]])) NA_real_ else tab[[x]][[cn]], numeric(1))
  cat(sprintf("PCA metrics found for %d of %d networks\n", sum(is.finite(NET$PC1_var)), nrow(NET)))
} else cat("\nNo *_PCA_variance.csv in the topo folder: PCA metrics skipped (they need RACIPE runs).\n")

write.csv(NET,  file.path(outdir, "perturbation_network_level.csv"), row.names = FALSE)
write.csv(MODL, file.path(outdir, "perturbation_module_level.csv"),  row.names = FALSE)
write.csv(PAIR, file.path(outdir, "perturbation_pair_level.csv"),    row.names = FALSE)

## ---- validation against the existing Metrics tables -----------------------------------------------------
ONET <- "Comb_AS_C_E_M_RS"; o <- NET[NET$kind == "orig", ]
vf <- "Metrics/metrics_network_level.csv"
if (file.exists(vf)) {
  old <- read.csv(vf, stringsAsFactors = FALSE); old <- old[old$network == ONET & old$l_max == LMAX, ]
  if (nrow(old)) {
    cmp <- intersect(names(o), names(old)); cmp <- setdiff(cmp, c("network", "kind", "rep", "n_removed", "n_added", "n_states"))
    v <- data.frame(metric = cmp, this_run = sapply(cmp, function(c) o[[c]]), metrics_csv = sapply(cmp, function(c) old[[c]]))
    v$diff <- v$this_run - v$metrics_csv
    cat("\nVALIDATION of the original network against Metrics/metrics_network_level.csv (l_max", LMAX, "):\n")
    print(v, row.names = FALSE); write.csv(v, file.path(outdir, "validation_original.csv"), row.names = FALSE)
    cat("T_* , density and impurity should agree to rounding; P_s columns can differ because the\nBoolean run was repeated (basin sizes drift).\n")
  } else cat("\n(", ONET, "is not in", vf, ": no validation)\n")
} else cat("\n(no Metrics/metrics_network_level.csv found: no validation)\n")

## ---- figures ----------------------------------------------------------------------------------------------
KC <- c(orig = "black", lower = "#2a78d6", higher = "#e34948", dmatch = "#1baf7a")
KL <- c(orig = "original", lower = "remove edges", higher = "add edges", dmatch = "density-matched")
STRUCT <- c(T_network = "network team strength", T_inter_blocks_mean = "inter-module strength (team-resolved, mean)",
            T_inter_mean = "inter-module strength (|block mean|, mean)", T_intra_mean = "intra-module strength (mean)",
            density = "density", impurity = "impurity")
OUTC <- c(ps_net_infl_top30 = "P_s network, influence-weighted", ps_net_state_top30 = "P_s network, state-based",
          ps_mod_infl_top30 = "P_s module mean, influence-weighted", ps_mod_state_abs_top30 = "|P_s| module mean, state-based",
          PC1_var = "PC1 variance", NumPC_90 = "PCs to 90%", PCA_entropy_norm = "PCA entropy (norm)")
OUTC <- OUTC[names(OUTC) %in% names(NET) & sapply(names(OUTC), function(c) any(is.finite(NET[[c]])))]
STRUCT <- STRUCT[names(STRUCT) %in% names(NET)]
NET$kind_lab <- factor(KL[NET$kind], levels = unname(KL))
NEW <- NET[NET$kind != "orig", ]; ORIG <- NET[NET$kind == "orig", ]
thm <- function() theme_minimal(base_size = 10) + theme(plot.title = element_text(face = "bold"),
  plot.subtitle = element_text(colour = "#52514e", size = 8.5), plot.caption = element_text(colour = "#52514e", size = 7, hjust = 0),
  legend.position = "top", panel.grid.minor = element_blank(), strip.text = element_text(size = 7.5),
  plot.background = element_rect(fill = "white", colour = NA))
long <- function(d) do.call(rbind, lapply(names(OUTC), function(y) do.call(rbind, lapply(names(STRUCT), function(x)
  data.frame(network = d$network, kind_lab = d$kind_lab, x = d[[x]], y = d[[y]],
             xl = factor(STRUCT[[x]], levels = unname(STRUCT)), yl = factor(OUTC[[y]], levels = unname(OUTC)))))))
Ln <- long(NEW); Lo <- long(ORIG)
rl <- do.call(rbind, lapply(split(Ln, list(Ln$xl, Ln$yl)), function(d) if (nrow(d) > 2 && sd(d$x) > 0 && sd(d$y) > 0)
  data.frame(xl = d$xl[1], yl = d$yl[1], lab = sprintf("r = %+.2f", cor(d$x, d$y)))))
S1 <- ggplot(Ln, aes(x, y)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey30", linewidth = 0.5, formula = y ~ x) +
  geom_point(aes(colour = kind_lab), size = 1.6, alpha = 0.8) +
  geom_point(data = Lo, shape = 18, size = 3.2, colour = "black") +
  geom_text(data = rl, aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.4, size = 2.5, fontface = "bold", inherit.aes = FALSE) +
  facet_grid(yl ~ xl, scales = "free", labeller = labeller(.rows = label_wrap_gen(18), .cols = label_wrap_gen(18))) +
  scale_colour_manual(values = setNames(unname(KC[-1]), unname(KL[-1])), name = NULL) +
  labs(title = sprintf("Landscape against structure across the perturbed networks of %s", ONET),
       subtitle = "one point per network (30), black diamond = original; line = least-squares fit over the 30", x = NULL, y = NULL,
       caption = "P_s is the top-30% aggregate. The influence-weighted score shares its influence matrix with T_s (partly structural); the state-based score does not.") + thm()
ggsave(file.path(outdir, "scatter_landscape_vs_structure.png"), S1, width = 3 + 2.4 * length(STRUCT), height = 1.6 + 1.9 * length(OUTC), dpi = 250, bg = "white")

## correlation matrix: all 30, and each kind on its own (10 each: noisy)
subsets <- list(all = NEW, lower = NEW[NEW$kind == "lower", ], higher = NEW[NEW$kind == "higher", ], dmatch = NEW[NEW$kind == "dmatch", ])
cm <- do.call(rbind, lapply(names(subsets), function(s) do.call(rbind, lapply(names(OUTC), function(y) do.call(rbind, lapply(names(STRUCT), function(x) {
  d <- subsets[[s]]; ok <- is.finite(d[[x]]) & is.finite(d[[y]])
  r <- if (sum(ok) > 2 && sd(d[[x]][ok]) > 0 && sd(d[[y]][ok]) > 0) cor(d[[x]][ok], d[[y]][ok]) else NA_real_
  data.frame(subset = sprintf("%s (n = %d)", c(all = "all perturbed", lower = "remove edges", higher = "add edges", dmatch = "density-matched")[[s]], nrow(d)),
             x = factor(STRUCT[[x]], levels = unname(STRUCT)), y = factor(OUTC[[y]], levels = rev(unname(OUTC))), r = r) }))))))
cm$subset <- factor(cm$subset, levels = unique(cm$subset))
CM <- ggplot(cm, aes(x, y, fill = r)) + geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = ifelse(is.finite(r), sprintf("%.2f", r), "n/a"), colour = is.finite(r) & abs(r) > 0.6), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
  scale_fill_gradient2(low = "#2a78d6", mid = "#f0efec", high = "#e34948", midpoint = 0, limits = c(-1, 1), na.value = "grey95", name = "Pearson r") +
  facet_wrap(~subset, nrow = 1) +
  labs(title = "Correlation of landscape with structure across the perturbed networks",
       subtitle = "left: all 30; right: each kind on its own (10 networks: read the pattern, not single values)", x = NULL, y = NULL,
       caption = "Within a kind the structural range is small and basin sizes carry simulation noise, so single correlations there are weak evidence.") +
  thm() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1, size = 7.5), axis.text.y = element_text(size = 8))
ggsave(file.path(outdir, "correlation_matrix.png"), CM, width = 17, height = 1.6 + 0.65 * length(OUTC), dpi = 250, bg = "white")

## the three kinds against the original
KV <- do.call(rbind, lapply(names(OUTC), function(y) data.frame(kind_lab = NET$kind_lab, v = NET[[y]], yl = factor(OUTC[[y]], levels = unname(OUTC)))))
ref <- KV[KV$kind_lab == KL[["orig"]], ]
nK <- length(OUTC); colsK <- ceiling(nK / 2)
K1 <- ggplot(KV, aes(kind_lab, v, colour = kind_lab)) +
  geom_hline(data = ref, aes(yintercept = v), linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_boxplot(data = KV[KV$kind_lab != KL[["orig"]], ], aes(fill = kind_lab), alpha = 0.18, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.8) +
  facet_wrap(~yl, scales = "free_y", nrow = 2, labeller = label_wrap_gen(26)) +
  scale_x_discrete(limits = unname(KL), labels = c("original", "remove\nedges", "add\nedges", "density-\nmatched")) +
  scale_colour_manual(values = setNames(unname(KC), unname(KL)), guide = "none") + scale_fill_manual(values = setNames(unname(KC), unname(KL)), guide = "none") +
  labs(title = "Landscape metrics of the perturbed networks against the original", subtitle = "each point is one network; dashed = original", x = NULL, y = NULL) +
  thm() + theme(axis.text.x = element_text(size = 8), strip.text = element_text(size = 8.5), panel.spacing = unit(1.1, "lines"))
ggsave(file.path(outdir, "kinds_vs_original.png"), K1, width = 3.6 * colsK + 1, height = 7.2, dpi = 250, bg = "white")

cat("\nmean (sd) by kind:\n")
for (c in c(names(STRUCT), names(OUTC))) {
  a <- tapply(NET[[c]], NET$kind, function(x) sprintf("%.3f (%.3f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE)))
  cat(sprintf("  %-26s orig %s | lower %s | higher %s | dmatch %s\n", c, sprintf("%.3f", ORIG[[c]]), a[["lower"]], a[["higher"]], a[["dmatch"]]))
}
cat(sprintf("\nfiles written to %s\n", outdir))