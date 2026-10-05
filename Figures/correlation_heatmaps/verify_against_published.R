suppressMessages({library(dplyr)})
repo <- getwd()
sfx <- c("2"="", "3"="_triple", "4"="_four", "5"="_five")
dirs <- c("2"="Double_EMT_MecnSens/correlations_double%s",
          "3"="Triple_EMT_MechSens/correlations_triple%s",
          "4"="Four_EMT_MechSens/correlations_four%s",
          "5"="Five_EMT_MechSens/correlations_five%s")
pubfile <- function(k, v) {
  d <- sprintf(dirs[[k]], if (v=="state") "_stateps" else "")
  file.path(d, sprintf("correlations_pearson_equal_lmax10%s.csv", sfx[[k]]))
}
mine <- list(infl  = read.csv("Figures/correlation_heatmaps/corr_heatmap_infl_lmax10.csv"),
             state = read.csv("Figures/correlation_heatmaps/corr_heatmap_state_lmax10.csv"))
# published row -> my output label ; published col -> my input name
ROW <- c(dominant_module="module P_s (dominant)", top30_module="module P_s (top30)",
         top60_module="module P_s (top60)", dominant_network="network P_s (dominant)",
         top30_network="network P_s (top30)", top60_network="network P_s (top60)",
         PC1_var="PC1 variance", NumPC_90="PCs to 90%", PCA_entropy_norm="PCA entropy (norm)")
COL <- c(network_T_s="T_network", intra_module_sum_T_s="T_network_intra_sum",
         inter_module_sum_T_s="T_inter_sum", impurity="impurity", density="density",
         intra_density="intra_density", inter_density="inter_density")
all <- list()
for (v in names(mine)) for (k in names(dirs)) {
  f <- pubfile(k, v); if (!file.exists(f)) { cat("missing:", f, "\n"); next }
  pub <- read.csv(f, stringsAsFactors = FALSE)
  m <- mine[[v]] %>% filter(n_modules == as.integer(k))
  for (r in intersect(names(ROW), pub$state_type)) for (cc in intersect(names(COL), names(pub))) {
    pv <- pub[pub$state_type == r, cc]
    mv <- m$r[m$output == ROW[[r]] & m$input == COL[[cc]]]
    if (!length(mv)) next
    all[[length(all)+1]] <- data.frame(variant=v, n=as.integer(k), row=r, col=cc,
                                       published=pv, mine=mv, diff=abs(pv-mv))
  }
}
A <- bind_rows(all)
cat(sprintf("\ncompared %d cells\n", nrow(A)))
cat(sprintf("match to <0.001 : %d\nmatch to <0.01  : %d\nmatch to <0.05  : %d\n",
    sum(A$diff<0.001), sum(A$diff<0.01), sum(A$diff<0.05)))
cat("\n--- by variant x n: median |diff| and worst ---\n")
print(A %>% group_by(variant,n) %>% summarise(cells=n(), med=round(median(diff),4),
      worst=round(max(diff),4), .groups="drop"), n=20)
cat("\n--- the 12 largest disagreements ---\n")
print(as.data.frame(A %>% arrange(desc(diff)) %>% head(12)), row.names=FALSE, digits=4)
