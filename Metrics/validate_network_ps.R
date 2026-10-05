## Validate the network-level INFLUENCE-WEIGHTED P_s in metrics_network_level.csv against the
## dominant_network / top30_network / top60_network rows of the published
## correlations_pearson_equal_*.csv. Those coefficients are all that exists --
## the per-network values themselves were never saved.
repo <- getwd()
N <- read.csv(file.path(repo,"Metrics","metrics_network_level.csv"), stringsAsFactors=FALSE)
SRC <- list(c("Double_EMT_MecnSens","correlations_double",2),
            c("Triple_EMT_MechSens","correlations_triple",3),
            c("Four_EMT_MechSens","correlations_four",4),
            c("Five_EMT_MechSens","correlations_five",5))
## published column name -> column name in metrics_network_level.csv
NETCOL <- c(network_T_s = "T_network", impurity = "impurity", density = "density",
            intra_density = "intra_density", inter_density = "inter_density")
ok <- c(dominant=0, top30=0, top60=0); tot <- ok; worst <- ok
for (S in SRC) for (lm in c(1,5,7,10)) {
  cf <- Sys.glob(file.path(repo,S[1],S[2],"correlations_pearson_equal_*.csv"))
  cf <- cf[grepl(paste0("lmax",lm,"(_|\\.)"), cf)]; if (length(cf)!=1) next
  pub <- read.csv(cf, stringsAsFactors=FALSE); rownames(pub) <- pub$state_type
  D <- N[N$n_modules == as.integer(S[3]) & N$l_max == lm, ]
  if (!nrow(D)) next
  for (v in names(ok)) for (pc in names(NETCOL)) {
    want <- pub[paste0(v,"_network"), pc]; if (is.null(want) || is.na(want)) next
    got  <- cor(D[[paste0("P_s_net_infl_", v)]], D[[NETCOL[[pc]]]])
    tot[v] <- tot[v] + 1
    d <- abs(got - want); worst[v] <- max(worst[v], d)
    if (d < 6e-5) ok[v] <- ok[v] + 1
  }
}
cat("\n=== network-level P_s vs published *_network correlation rows\n")
for (v in names(ok))
  cat(sprintf("  %-9s %3d / %3d reproduced   worst deviation %.2e\n", v, ok[v], tot[v], worst[v]))
