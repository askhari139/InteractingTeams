## compute_inter_teamaware.R ----------------------------------------------
## Team-aware inter-module strength, written to a SEPARATE file so the
## validated metrics_*.csv lineage (T_inter reproduces the published values
## 1800/1800) is left untouched.
##
##   Metrics/metrics_inter_teamaware.csv   network x ordered pair x l_max
##   Metrics/metrics_inter_teamaware_network.csv   network x l_max roll-ups
##
## Usage:  Rscript Metrics/compute_inter_teamaware.R [lmax1,lmax2,...]
## -------------------------------------------------------------------------
suppressMessages(library(dplyr))
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))
args <- commandArgs(trailingOnly = TRUE)
LMAX <- if (length(args)) as.integer(strsplit(args[1], ",")[[1]]) else c(1, 5, 7, 10)

N <- read.csv(file.path(repo, "Metrics", "metrics_network_level.csv"))
nets <- sort(unique(N$network))
cat("networks:", length(nets), " l_max:", paste(LMAX, collapse = ","), "\n")

rows <- list()
for (net in nets) {
  mods <- modules_of(net); tp <- find_net_file(repo, net, "topo")
  mt <- setNames(lapply(mods, function(m) module_teams(repo, m)), mods)
  mn <- setNames(lapply(mods, function(m) module_nodes(repo, m)), mods)
  for (lm in LMAX) {
    old <- setwd(dirname(tp)); I <- influence_full(basename(tp), lm); setwd(old)
    for (a in mods) for (b in mods) if (a != b) {
      rows[[length(rows) + 1]] <- data.frame(
        network = net, l_max = lm, n_modules = length(mods),
        module_from = a, module_to = b,
        T_inter_old     = T_inter(I, mn[[a]], mn[[b]]),
        T_inter_blocks  = T_inter_blocks(I, mt[[a]], mt[[b]]),
        T_inter_aligned = T_inter_aligned(I, mt[[a]], mt[[b]]),
        stringsAsFactors = FALSE)
    }
  }
  cat(".", sep = "")
}
cat("\n")
P <- bind_rows(rows)
write.csv(P, file.path(repo, "Metrics", "metrics_inter_teamaware.csv"), row.names = FALSE)

NW <- P %>% group_by(network, l_max, n_modules) %>%
  summarise(T_inter_sum_old          = sum(T_inter_old),
            T_inter_blocks_sum       = sum(T_inter_blocks),
            T_inter_blocks_mean      = mean(T_inter_blocks),
            T_inter_aligned_sum      = sum(T_inter_aligned),
            T_inter_aligned_abs_sum  = sum(abs(T_inter_aligned)),
            T_inter_aligned_mean     = mean(T_inter_aligned),
            ## coherence: how much of the block magnitude is team-aligned
            T_inter_coherence        = sum(abs(T_inter_aligned)) / sum(T_inter_blocks),
            .groups = "drop")
write.csv(NW, file.path(repo, "Metrics", "metrics_inter_teamaware_network.csv"), row.names = FALSE)
cat(sprintf("wrote Metrics/metrics_inter_teamaware.csv          %d rows\n", nrow(P)))
cat(sprintf("wrote Metrics/metrics_inter_teamaware_network.csv  %d rows\n", nrow(NW)))
