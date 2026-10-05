## validate_metrics.R ------------------------------------------------------
## Recompute every metric from the .topo / .teams / finFlagFreq inputs and
## compare against Combined_data/combined_module_level_data.csv, which is the
## compiled form of the existing correlations_* outputs.
## -------------------------------------------------------------------------

repo <- getwd()
options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

truth <- read.csv(file.path(repo, "Combined_data", "combined_module_level_data.csv"),
                  stringsAsFactors = FALSE)
truth <- truth[truth$pheno_score_method == "newNorm", ]

nets  <- sort(unique(truth$network_name))
lmaxs <- sort(unique(truth$l_max))

tally <- new.env(parent = emptyenv())
fails <- list()
note  <- function(metric, ok, label = NA, got = NA, want = NA) {
  ok <- isTRUE(ok)
  cur <- mget(metric, tally, ifnotfound = list(c(0, 0)))[[1]]
  assign(metric, cur + c(ok, 1), envir = tally)
  if (!ok) fails[[length(fails) + 1]] <<-
    sprintf("%-12s %-28s got=%s want=%s", metric, label, format(got), format(want))
}
TOL <- 1e-9

for (net in nets) {
  topo_p  <- find_net_file(repo, net, "topo")
  teams_p <- find_net_file(repo, net, "teams")
  topo    <- read_topo(topo_p)
  teams   <- read_teams(teams_p)
  mods    <- modules_of(net)
  modn    <- setNames(lapply(mods, function(m) module_nodes(repo, m)), mods)
  modt    <- setNames(lapply(mods, function(m) module_teams(repo, m)), mods)

  ## topology-only metrics (lmax independent)
  ref <- truth[truth$network_name == net, ][1, ]
  ds  <- density_split(topo, teams)
  note("density",       abs(density_core(topo, teams) - ref$density)       < TOL)
  note("impurity",      abs(impurity(topo, teams)     - ref$impurity)      < TOL)
  note("intra_density", abs(ds[["intra"]]             - ref$intra_density) < TOL)
  note("inter_density", abs(ds[["inter"]]             - ref$inter_density) < TOL)

  for (lm in lmaxs) {
    old <- setwd(dirname(topo_p))
    I   <- influence_full(basename(topo_p), lm)
    setwd(old)

    sub <- truth[truth$network_name == net & truth$l_max == lm, ]
    tn <- T_network(I, teams)
    note("network_T_s", abs(tn - sub$network_T_s[1]) < TOL,
         sprintf("%s lmax%d", net, lm), tn, sub$network_T_s[1])

    for (i in seq_len(nrow(sub))) {
      m <- sub$module_abbr[i]
      tm <- T_module_intra(I, modn[[m]], modt[[m]])
      note("module_T_s", abs(tm - sub$module_T_s[i]) < TOL,
           sprintf("%s/%s lmax%d", net, m, lm), tm, sub$module_T_s[i])
    }
  }
}

cat("\n=== recomputed vs existing outputs (", length(nets), "networks x",
    length(lmaxs), "lmax )\n")
for (k in sort(ls(tally))) {
  v <- get(k, tally)
  cat(sprintf("  %-14s %4d / %4d match%s\n", k, v[1], v[2],
              ifelse(v[1] == v[2], "", "   <-- MISMATCH")))
}
if (length(fails)) {
  cat("\nfirst mismatches:\n")
  for (f in head(unlist(fails), 15)) cat("  ", f, "\n")
}
