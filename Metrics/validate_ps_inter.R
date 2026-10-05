## Validate (a) inter-module team strength against the panel files and
## (b) the phenotypic score against the module_data P_s columns.
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

PANEL <- list(
  c("Double_EMT_MecnSens","correlations_double"),
  c("Triple_EMT_MechSens","correlations_triple"),
  c("Four_EMT_MechSens","correlations_four"),
  c("Five_EMT_MechSens","correlations_five"))

## ---- (a) inter-module team strength --------------------------------------
ok <- 0; tot <- 0; bad <- character()
for (p in PANEL) {
  d <- file.path(repo, p[1], p[2], "Team_strength_pearson_equal")
  for (f in list.files(d, "panel", full.names = TRUE)) {
    lm <- as.integer(sub(".*lmax(\\d+)_.*", "\\1", basename(f)))
    pan <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(pan))) {
      net <- pan$network_name[i]
      tp  <- find_net_file(repo, net, "topo")
      old <- setwd(dirname(tp)); I <- influence_full(basename(tp), lm); setwd(old)
      mods <- modules_of(net)
      for (a in mods) for (b in mods) if (a != b) {
        col <- paste0("T_", a, "→", b)
        if (!col %in% names(pan)) next
        want <- pan[[col]][i]
        if (is.na(want)) next
        got <- T_inter(I, module_nodes(repo, a), module_nodes(repo, b))
        tot <- tot + 1
        if (!is.na(got) && abs(got - want) < 1e-9) ok <- ok + 1
        else bad <- c(bad, sprintf("%s %s lmax%d got=%.10f want=%.10f", net, col, lm, got, want))
      }
    }
  }
}
cat(sprintf("\ninter-module T_A->B: %d / %d match\n", ok, tot))
for (b in head(bad, 6)) cat("   ", b, "\n")

## ---- (b) phenotypic score -------------------------------------------------
truth <- read.csv(file.path(repo,"Combined_data","combined_module_level_data.csv"),
                  stringsAsFactors = FALSE)
freq_file <- function(net) {
  for (d in c("Two_module_nets", NETWORK_DIR)) {
    p <- file.path(repo, d, paste0(net, "_finFlagFreq.csv"))
    if (file.exists(p)) return(p)
  }
  NA_character_
}
cat("\nphenotypic score, first few networks (raw formula-sheet dF):\n")
cat(sprintf("  %-16s %-4s %10s %10s %12s %12s\n","network","mod","got_dom","newNorm","stateps","ratio(sp/got)"))
done <- 0
for (net in unique(truth$network_name)) {
  fp <- freq_file(net); if (is.na(fp)) next
  tp <- find_net_file(repo, net, "topo")
  topo <- read_topo(tp)
  nds  <- read_nodes(file.path(dirname(tp), paste0(net, "_nodes.txt")))
  for (m in modules_of(net)) {
    mt <- restrict_teams(module_teams(repo, m), nds)
    if (length(mt) < 2) next
    df <- state_pheno_scores(fp, topo, nds, mt)
    agg <- aggregate_ps(df)
    r1 <- truth[truth$network_name==net & truth$module_abbr==m &
                truth$l_max==1 & truth$pheno_score_method=="newNorm", ]
    r2 <- truth[truth$network_name==net & truth$module_abbr==m &
                truth$l_max==1 & truth$pheno_score_method=="stateps", ]
    cat(sprintf("  %-16s %-4s %10.5f %10.5f %12.5f %12.4f\n", net, m,
        agg[["dominant"]], r1$P_s_dominant[1], r2$P_s_dominant[1],
        r2$P_s_dominant[1]/agg[["dominant"]]))
  }
  done <- done + 1; if (done >= 4) break
}
