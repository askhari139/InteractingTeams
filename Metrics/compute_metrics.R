## compute_metrics.R -------------------------------------------------------
## Compute every metric on the formula sheet for all Double/Triple/Four/Five
## networks and write two tidy CSVs:
##
##   Metrics/metrics_network_level.csv  network x l_max
##   Metrics/metrics_module_level.csv   network x module x l_max
##   Metrics/metrics_pair_level.csv     network x (module A -> module B) x l_max
##
## Usage:  Rscript Metrics/compute_metrics.R  [lmax1,lmax2,...] [--singles]
##
## --singles computes the SINGLE-module networks instead: the Comb_<M>_<M>
## self-pairs in Two_module_nets, which are byte-identical in topology and
## node order to Single_<M>. They are normally excluded (see all_networks),
## which is why n_modules starts at 2 in the main tables. With the flag the
## output goes to Metrics/metrics_*_singles.csv so the validated 56-network
## tables are never touched.
## -------------------------------------------------------------------------

repo <- getwd()
options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

args    <- commandArgs(trailingOnly = TRUE)
SINGLES <- "--singles" %in% args
args    <- setdiff(args, "--singles")
LMAX    <- if (length(args)) as.integer(strsplit(args[1], ",")[[1]]) else c(1, 5, 7, 10)
SUFFIX  <- if (SINGLES) "_singles" else ""

all_networks <- function(repo) {
  nets <- unlist(lapply(NETWORK_DIR, function(d)
    sub("\\.topo$", "", basename(Sys.glob(file.path(repo, d, "Comb_*.topo"))))))
  ## drop self-pairs (single modules) and anything with an unknown module
  known <- sub("^Single_", "", sub("_nodes\\.txt$", "",
            basename(Sys.glob(file.path(repo, MODULE_DIR, "Single_*_nodes.txt")))))
  ## keep only networks whose every module has a genuine two-team split;
  ## A, CC, GBM, GP, Gc, Gm, OL and Tbs have a single NewProgram in
  ## allNodes.nodes, so no team strength or phenotypic score is defined for them
  nets <- nets[vapply(nets, function(n) {
    m <- modules_of(n)
    (if (SINGLES) length(m) == 1 else length(m) > 1) &&
      all(m %in% known) && all(vapply(m, function(x) has_two_teams(repo, x), logical(1)))
  }, logical(1))]
  sort(unique(nets))
}

freq_file <- function(repo, net) {
  for (d in c("Two_module_nets", NETWORK_DIR)) {
    p <- file.path(repo, d, paste0(net, "_finFlagFreq.csv"))
    if (file.exists(p)) return(p)
  }
  NA_character_
}

net_rows  <- list()
mod_rows  <- list()
pair_rows <- list()

nets <- all_networks(repo)
cat("networks:", length(nets), " l_max:", paste(LMAX, collapse = ","), "\n")

for (net in nets) {
  topo_p  <- find_net_file(repo, net, "topo")
  teams_p <- find_net_file(repo, net, "teams")
  topo    <- read_topo(topo_p)
  teams   <- read_teams(teams_p)
  mods    <- modules_of(net)
  modn    <- setNames(lapply(mods, function(m) module_nodes(repo, m)), mods)
  modt    <- setNames(lapply(mods, function(m) module_teams(repo, m)), mods)

  core <- core_nodes(topo)
  pc   <- pca_metrics(repo, net)
  dens <- density_split(topo, teams)
  topo_metrics <- list(impurity      = impurity(topo, teams),
                       density       = density_core(topo, teams),
                       intra_density = dens[["intra"]],
                       inter_density = dens[["inter"]])

  ## state-based score: depends only on the states, so computed once.
  ## Signed, basin-weighted *mean* over the top 30%/60% of the state space.
  fp   <- freq_file(repo, net)
  nds  <- read_nodes(file.path(dirname(topo_p), paste0(net, "_nodes.txt")))
  ps_state <- list()
  if (!is.na(fp)) {
    for (m in mods) {
      mt <- restrict_teams(modt[[m]], nds)
      if (length(mt) < 2) next  # no two-team structure
      ps_state[[m]] <- aggregate_ps(state_ps_table(fp, nds, mt), normalise = TRUE)
    }
  }

  for (lm in LMAX) {
    old <- setwd(dirname(topo_p)); I <- influence_full(basename(topo_p), lm); setwd(old)

    T_net   <- T_network(I, teams)
    T_intra <- setNames(lapply(mods, function(m)
                 T_module_intra(I, modn[[m]], modt[[m]])), mods)

    ## influence-weighted score: l_max dependent, so recomputed per l_max,
    ## as |dF| with renormalised basin weights (reproduces the published column)
    ps_infl <- list()
    if (!is.na(fp)) {
      for (m in mods) {
        mt <- restrict_teams(modt[[m]], nds)
        if (length(mt) < 2) next
        ps_infl[[m]] <- aggregate_ps(state_pheno_scores(fp, I, nds, mt),
                                     absolute = TRUE, normalise = TRUE)
      }
    }

    ## network-level phenotypic scores, network's own teams.
    ## influence-weighted: |dF| on the core matrix, l_max dependent.
    ## state-based: signed team-mean difference, l_max free.
    ps_net <- NULL; ps_nets <- NULL
    if (!is.na(fp)) {
      tb <- network_ps_table(fp, I, nds, teams, core)
      if (!is.null(tb)) ps_net <- aggregate_ps(tb, absolute = TRUE, normalise = TRUE)
      tbs <- network_state_ps_table(fp, nds, teams)
      if (!is.null(tbs)) ps_nets <- aggregate_ps(tbs, normalise = TRUE)
    }
    net_rows[[length(net_rows) + 1]] <- data.frame(
      network = net, l_max = lm, n_modules = length(mods),
      has_PS = as.integer("PS" %in% mods),
      n_nodes = length(nds), n_core = length(core),
      T_network = T_net,
      T_network_intra_sum = sum(unlist(T_intra), na.rm = TRUE),
      T_inter_sum  = sum(vapply(mods, function(a) sum(vapply(setdiff(mods, a),
                       function(b) T_inter(I, modn[[a]], modn[[b]]), numeric(1))),
                       numeric(1))),
      T_inter_mean = mean(unlist(lapply(mods, function(a) vapply(setdiff(mods, a),
                       function(b) T_inter(I, modn[[a]], modn[[b]]), numeric(1))))),
      PC1_var          = pc$PC1_var,
      NumPC_90         = pc$NumPC_90,
      PCA_entropy      = pc$PCA_entropy,
      PCA_entropy_norm = pc$PCA_entropy_norm,
      impurity = topo_metrics$impurity, density = topo_metrics$density,
      intra_density = topo_metrics$intra_density,
      inter_density = topo_metrics$inter_density,
      P_s_net_infl_dominant  = if (!is.null(ps_net))  ps_net[["dominant"]]  else NA_real_,
      P_s_net_infl_top30     = if (!is.null(ps_net))  ps_net[["top30"]]     else NA_real_,
      P_s_net_infl_top60     = if (!is.null(ps_net))  ps_net[["top60"]]     else NA_real_,
      P_s_net_state_dominant = if (!is.null(ps_nets)) ps_nets[["dominant"]] else NA_real_,
      P_s_net_state_top30    = if (!is.null(ps_nets)) ps_nets[["top30"]]    else NA_real_,
      P_s_net_state_top60    = if (!is.null(ps_nets)) ps_nets[["top60"]]    else NA_real_,
      stringsAsFactors = FALSE)

    for (m in mods) {
      others <- setdiff(mods, m)
      mod_rows[[length(mod_rows) + 1]] <- data.frame(
        network = net, l_max = lm, n_modules = length(mods), module = m,
        has_PS = as.integer("PS" %in% mods),
        T_network = T_net,
        T_module_intra = T_intra[[m]],
        T_module_out_sum = sum(vapply(others, function(o)
          T_inter(I, modn[[m]], modn[[o]]), numeric(1)), na.rm = TRUE),
        T_module_in_sum  = sum(vapply(others, function(o)
          T_inter(I, modn[[o]], modn[[m]]), numeric(1)), na.rm = TRUE),
        impurity = topo_metrics$impurity, density = topo_metrics$density,
        intra_density = topo_metrics$intra_density,
        inter_density = topo_metrics$inter_density,
        P_s_state_dominant = if (!is.null(ps_state[[m]])) ps_state[[m]][["dominant"]] else NA_real_,
        P_s_state_top30    = if (!is.null(ps_state[[m]])) ps_state[[m]][["top30"]]    else NA_real_,
        P_s_state_top60    = if (!is.null(ps_state[[m]])) ps_state[[m]][["top60"]]    else NA_real_,
        P_s_infl_dominant  = if (!is.null(ps_infl[[m]]))  ps_infl[[m]][["dominant"]]  else NA_real_,
        P_s_infl_top30     = if (!is.null(ps_infl[[m]]))  ps_infl[[m]][["top30"]]     else NA_real_,
        P_s_infl_top60     = if (!is.null(ps_infl[[m]]))  ps_infl[[m]][["top60"]]     else NA_real_,
        stringsAsFactors = FALSE)
    }

    for (a in mods) for (b in mods) if (a != b) {
      tab_ab <- T_inter(I, modn[[a]], modn[[b]])
      tab_ba <- T_inter(I, modn[[b]], modn[[a]])
      pair_rows[[length(pair_rows) + 1]] <- data.frame(
        network = net, l_max = lm, n_modules = length(mods),
        module_from = a, module_to = b,
        T_from_to = tab_ab, T_to_from = tab_ba,
        T_inter_sum  = tab_ab + tab_ba,
        T_intra_sum  = T_intra[[a]] + T_intra[[b]],
        T_intra_mean = (T_intra[[a]] + T_intra[[b]]) / 2,
        stringsAsFactors = FALSE)
    }
  }
  cat(".", sep = "")
}
cat("\n")

netl <- do.call(rbind, net_rows)
mod  <- do.call(rbind, mod_rows)
pair <- do.call(rbind, pair_rows)
netl <- netl[order(netl$n_modules, netl$network, netl$l_max), ]
mod  <- mod[order(mod$n_modules, mod$network, mod$l_max, mod$module), ]
if (!is.null(pair))   ## single-module networks produce no ordered pairs
  pair <- pair[order(pair$n_modules, pair$network, pair$l_max,
                     pair$module_from, pair$module_to), ]

fn <- function(base) file.path(repo, "Metrics", sprintf("metrics_%s%s.csv", base, SUFFIX))
write.csv(netl, fn("network_level"), row.names = FALSE)
write.csv(mod,  fn("module_level"),  row.names = FALSE)
cat(sprintf("wrote %s %d rows x %d cols\n", fn("network_level"), nrow(netl), ncol(netl)))
cat(sprintf("wrote %s %d rows x %d cols\n", fn("module_level"),  nrow(mod),  ncol(mod)))
if (!is.null(pair)) {   ## single-module networks have no ordered pairs
  write.csv(pair, fn("pair_level"), row.names = FALSE)
  cat(sprintf("wrote %s %d rows x %d cols\n", fn("pair_level"), nrow(pair), ncol(pair)))
}
