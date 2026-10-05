## Validate the influence-weighted phenotypic score against the published
## (non-stateps) P_s columns, across all networks, modules and l_max.
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

truth <- read.csv(file.path(repo,"Combined_data","combined_module_level_data.csv"),
                  stringsAsFactors = FALSE)
truth <- truth[truth$pheno_score_method == "newNorm", ]
freq_file <- function(net) {
  for (d in c("Two_module_nets", NETWORK_DIR)) {
    p <- file.path(repo, d, paste0(net, "_finFlagFreq.csv")); if (file.exists(p)) return(p) }
  NA_character_
}
ok <- setNames(integer(3), c("dominant","top30","top60")); tot <- ok
flip <- ok; bad <- character()
for (net in unique(truth$network_name)) {
  fp <- freq_file(net); if (is.na(fp)) next
  tp  <- find_net_file(repo, net, "topo")
  nds <- read_nodes(file.path(dirname(tp), paste0(net, "_nodes.txt")))
  for (lm in sort(unique(truth$l_max))) {
    old <- setwd(dirname(tp)); I <- influence_full(basename(tp), lm); setwd(old)
    for (m in modules_of(net)) {
      if (!has_two_teams(repo, m)) next
      mt <- restrict_teams(module_teams(repo, m), nds)
      d <- state_pheno_scores(fp, I, nds, mt)
      d$ps <- abs(d$ps)                       # published column is |dF|
      dd <- d[order(-d$freq), ]; cum <- cumsum(dd$freq)
      wmean <- function(th) { k <- which(cum >= th)[1]; if (is.na(k)) k <- nrow(dd)
        sub <- dd[seq_len(k), ]; sum(sub$ps*sub$freq)/sum(sub$freq) }
      agg <- c(dominant = dd$ps[1], top30 = wmean(0.30), top60 = wmean(0.60))
      r <- truth[truth$network_name==net & truth$module_abbr==m & truth$l_max==lm, ]
      if (!nrow(r)) next
      for (k in names(ok)) {
        want <- r[[paste0("P_s_", k)]][1]; got <- agg[[k]]
        tot[k] <- tot[k] + 1
        if (abs(got - want) < 1e-9) ok[k] <- ok[k] + 1
        else if (abs(got + want) < 1e-9) flip[k] <- flip[k] + 1
        else if (length(bad) < 6)
          bad <- c(bad, sprintf("%-18s %-3s lmax%-3d %-8s got=%+.8f want=%+.8f", net,m,lm,k,got,want))
      }
    }
  }
}
cat("\n=== influence-weighted P_s as |dF|, basin-weighted MEAN over top X%\n")
for (k in names(ok))
  cat(sprintf("  P_s_%-9s %4d / %4d exact, %4d sign-flipped, %4d other\n",
      k, ok[k], tot[k], flip[k], tot[k]-ok[k]-flip[k]))
for (b in bad) cat("   ", b, "\n")
