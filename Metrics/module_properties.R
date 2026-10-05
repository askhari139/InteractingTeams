## module_properties.R ----------------------------------------------------
## Intrinsic properties of each two-team module, for regressing against the
## leave-one-module-out correlation impact.
##
## Writes Metrics/module_properties.csv (one row per module x l_max).
## -------------------------------------------------------------------------
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

MODS <- c("AS","C","E","M","PS","RS")
pole <- function(m) {                      # node -> 1 (team1) or 2 (team2)
  g <- module_teams(repo, m)
  setNames(rep(seq_along(g), lengths(g)), unlist(g))
}
POLE <- do.call(c, lapply(MODS, pole))

## ---- within-module, l_max independent ----------------------------------
within_topo <- list()
for (m in MODS) {
  tp   <- file.path(repo, MODULE_DIR, paste0("Single_", m, ".topo"))
  topo <- read_topo(tp)
  g    <- module_teams(repo, m)
  core <- core_nodes(topo)
  ## edges of the module's own core subgraph, self-loops dropped
  e <- topo[topo$Source %in% core & topo$Target %in% core &
            topo$Source != topo$Target, , drop = FALSE]
  lab  <- POLE[c(e$Source)]; lab2 <- POLE[c(e$Target)]
  same <- lab == lab2
  n    <- length(core); sz <- lengths(restrict_teams(g, core))
  within_topo[[m]] <- data.frame(
    module = m,
    n_nodes = length(module_nodes(repo, m)), n_core = n,
    density_within        = nrow(e) / (n * (n - 1)),
    impurity_within       = sum((same & e$Type < 0) | (!same & e$Type > 0)) / nrow(e),
    intra_density_within  = sum(same)  / sum(sz * (sz - 1)),
    inter_density_within  = sum(!same) / (2 * prod(sz)),
    stringsAsFactors = FALSE)
}
within_topo <- do.call(rbind, within_topo)

## ---- across-module, pooled over every network containing the module -----
known <- MODS
nets <- unlist(lapply(NETWORK_DIR, function(d)
  sub("\\.topo$", "", basename(Sys.glob(file.path(repo, d, "Comb_*.topo"))))))
nets <- sort(unique(nets[vapply(nets, function(n) {
  mm <- modules_of(n); length(mm) > 1 && all(mm %in% known) }, logical(1))]))

acc <- setNames(lapply(MODS, function(m) list(e = 0, mis = 0, pairs = 0, nets = 0)), MODS)
for (net in nets) {
  topo <- read_topo(find_net_file(repo, net, "topo"))
  mm   <- modules_of(net)
  nodesOf <- setNames(lapply(mm, function(x) module_nodes(repo, x)), mm)
  owner <- unlist(lapply(mm, function(x) setNames(rep(x, length(nodesOf[[x]])), nodesOf[[x]])))
  e <- topo[topo$Source %in% names(owner) & topo$Target %in% names(owner), , drop = FALSE]
  oa <- owner[c(e$Source)]; ob <- owner[c(e$Target)]
  cross <- oa != ob
  ec <- e[cross, , drop = FALSE]; oa <- oa[cross]; ob <- ob[cross]
  pa <- POLE[c(ec$Source)]; pb <- POLE[c(ec$Target)]
  ## "aligned": same pole & activating, or opposite pole & inhibiting
  mis <- !((pa == pb & ec$Type > 0) | (pa != pb & ec$Type < 0))
  for (m in mm) {
    sel <- oa == m | ob == m
    nm  <- length(nodesOf[[m]]); no <- sum(lengths(nodesOf)) - nm
    acc[[m]]$e     <- acc[[m]]$e     + sum(sel)
    acc[[m]]$mis   <- acc[[m]]$mis   + sum(mis[sel])
    acc[[m]]$pairs <- acc[[m]]$pairs + 2 * nm * no
    acc[[m]]$nets  <- acc[[m]]$nets  + 1
  }
}
across <- do.call(rbind, lapply(MODS, function(m) data.frame(
  module = m, n_networks = acc[[m]]$nets,
  cross_edges      = acc[[m]]$e,
  density_across   = acc[[m]]$e   / acc[[m]]$pairs,
  impurity_across  = if (acc[[m]]$e > 0) acc[[m]]$mis / acc[[m]]$e else NA_real_,
  stringsAsFactors = FALSE)))

## ---- intrinsic team strength, per l_max --------------------------------
ts <- list()
for (lm in c(1, 5, 7, 10)) {
  for (m in MODS) {
    old <- setwd(file.path(repo, MODULE_DIR))
    I <- influence_full(paste0("Single_", m, ".topo"), lm)
    setwd(old)
    ts[[length(ts) + 1]] <- data.frame(module = m, l_max = lm,
      T_s_module = T_module_intra(I, module_nodes(repo, m), module_teams(repo, m)),
      stringsAsFactors = FALSE)
  }
}
ts <- do.call(rbind, ts)

P <- merge(merge(ts, within_topo, by = "module"), across, by = "module")
P <- P[order(P$l_max, P$module), ]
write.csv(P, file.path(repo, "Metrics", "module_properties.csv"), row.names = FALSE)
cat(sprintf("wrote Metrics/module_properties.csv  %d rows x %d cols\n", nrow(P), ncol(P)))
print(P[P$l_max == 10, ], row.names = FALSE, digits = 4)
