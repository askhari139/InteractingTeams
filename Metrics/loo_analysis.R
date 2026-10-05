## loo_analysis.R ---------------------------------------------------------
## Leave-one-module-out effect on correlations, per analysis.md line 10.
##
## For every (outcome measure, predictor, l_max):
##   r_all       = cor over all 56 networks
##   r_without_X = cor over the 26 networks that do not contain module X
##   impact      = r_all - r_without_X
## plus a null: r over 2000 random 26-of-56 subsets, since a module removal
## always leaves 26 networks and r on 26 points is noisy.
##
## Outcome measures, in four families:
##   psM_infl   module-level influence-weighted P_s (top30, top60), averaged
##              over the network's modules -- the aggregation that reproduces
##              the published *_module correlation rows
##   psM_state  module-level state-based P_s (top30, top60), same aggregation
##   psN_infl   network-level influence-weighted P_s (top30, top60), the
##              network's own teams on the core matrix -- |dF|
##   psN_state  network-level state-based P_s (top30, top60), signed
##   pc         PCA metrics (PC1_var, NumPC_90, PCA_entropy, PCA_entropy_norm)
##
## Predictors: T_network and T_inter_sum (total inter-module team strength).
## -------------------------------------------------------------------------
repo <- getwd()
M    <- read.csv(file.path(repo,"Metrics","metrics_module_level.csv"),  stringsAsFactors=FALSE)
NET  <- read.csv(file.path(repo,"Metrics","metrics_network_level.csv"), stringsAsFactors=FALSE)

MODS <- c("AS","C","E","M","PS","RS")
PREDICTORS <- c(T_network = "network T_s", T_inter_sum = "inter-module T_s sum")

## module-level measures are averaged over the network's modules
MOD_MEASURE <- c(psM_infl_top30  = "P_s_infl_top30",  psM_infl_top60  = "P_s_infl_top60",
                 psM_state_top30 = "P_s_state_top30", psM_state_top60 = "P_s_state_top60")
## network-level measures are read straight off the network table
NET_MEASURE <- c(psN_infl_top30  = "P_s_net_infl_top30",
                 psN_infl_top60  = "P_s_net_infl_top60",
                 psN_state_top30 = "P_s_net_state_top30",
                 psN_state_top60 = "P_s_net_state_top60",
                 pc_PC1_var = "PC1_var", pc_NumPC_90 = "NumPC_90",
                 pc_PCA_entropy = "PCA_entropy", pc_PCA_entropy_norm = "PCA_entropy_norm")
FAMILY <- function(m) sub("_(top30|top60|PC1_var|NumPC_90|PCA_entropy|PCA_entropy_norm)$", "", m)

modsOf <- lapply(split(M$module, M$network), unique)

## one network-level table of every outcome and predictor, per l_max
build <- function(lm) {
  Ml <- M[M$l_max == lm, ]; Nl <- NET[NET$l_max == lm, ]
  d <- Nl[, c("network", names(PREDICTORS))]
  for (k in names(MOD_MEASURE))
    d[[k]] <- as.numeric(tapply(Ml[[MOD_MEASURE[[k]]]], Ml$network, mean)[d$network])
  for (k in names(NET_MEASURE)) d[[k]] <- Nl[[NET_MEASURE[[k]]]][match(d$network, Nl$network)]
  d
}
MEASURES <- c(names(MOD_MEASURE), names(NET_MEASURE))

rows <- list(); nl <- list()
set.seed(1); NDRAW <- 2000
for (lm in sort(unique(M$l_max))) {
  d <- build(lm)
  for (pv in names(PREDICTORS)) for (mv in MEASURES) {
    x <- d[[pv]]; y <- d[[mv]]
    keepAll <- is.finite(x) & is.finite(y)
    r_all <- cor(y[keepAll], x[keepAll])
    rows[[length(rows)+1]] <- data.frame(l_max=lm, predictor=pv, measure=mv,
        family=FAMILY(mv), removed="All", n=sum(keepAll), r=r_all, impact=0,
        stringsAsFactors=FALSE)
    for (X in MODS) {
      k <- keepAll & !vapply(d$network, function(n) X %in% modsOf[[n]], logical(1))
      r <- cor(y[k], x[k])
      rows[[length(rows)+1]] <- data.frame(l_max=lm, predictor=pv, measure=mv,
          family=FAMILY(mv), removed=X, n=sum(k), r=r, impact=r_all-r,
          stringsAsFactors=FALSE)
    }
    idx <- which(keepAll)
    rr <- replicate(NDRAW, { i <- sample(idx, 26); cor(y[i], x[i]) })
    nl[[length(nl)+1]] <- data.frame(l_max=lm, predictor=pv, measure=mv,
        family=FAMILY(mv), draw=seq_len(NDRAW), r=rr, r_all=r_all,
        stringsAsFactors=FALSE)
  }
}
L <- do.call(rbind, rows); NULLD <- do.call(rbind, nl)
write.csv(L,     file.path(repo,"Metrics","loo_correlations.csv"), row.names=FALSE)
write.csv(NULLD, file.path(repo,"Metrics","loo_null.csv"),         row.names=FALSE)

LX <- L[L$removed != "All", ]

cat("\n=== baseline r over all 56 networks\n")
b <- L[L$removed=="All", c("predictor","measure","l_max","r")]
print(reshape(b, idvar=c("predictor","measure"), timevar="l_max", direction="wide"),
      row.names=FALSE, digits=3)

cat("\n=== mean |impact| by module, per predictor and family\n")
agg <- aggregate(abs(LX$impact),
                 by=list(predictor=LX$predictor, family=LX$family, module=LX$removed), FUN=mean)
names(agg)[4] <- "mean_abs_impact"
for (pv in names(PREDICTORS)) for (fm in unique(LX$family)) {
  s <- agg[agg$predictor==pv & agg$family==fm, c("module","mean_abs_impact")]
  s <- s[order(-s$mean_abs_impact), ]
  cat(sprintf("  %-12s %-10s : %s\n", pv, fm,
      paste(sprintf("%s=%.3f", s$module, s$mean_abs_impact), collapse="  ")))
}

cat("\n=== observed |impact| vs the random-subset null\n")
cmp <- list()
for (lm in sort(unique(L$l_max))) for (pv in names(PREDICTORS)) for (mv in MEASURES) {
  nd  <- NULLD[NULLD$l_max==lm & NULLD$predictor==pv & NULLD$measure==mv, ]
  p95 <- as.numeric(quantile(abs(nd$r_all[1] - nd$r), .95))
  obs <- LX[LX$l_max==lm & LX$predictor==pv & LX$measure==mv, ]
  cmp[[length(cmp)+1]] <- data.frame(l_max=lm, predictor=pv, measure=mv,
    null_p95=p95, obs_max=max(abs(obs$impact)), n_above=sum(abs(obs$impact) > p95))
}
CMP <- do.call(rbind, cmp)
write.csv(CMP, file.path(repo,"Metrics","loo_null_summary.csv"), row.names=FALSE)
sm <- aggregate(cbind(null_p95, obs_max, n_above) ~ predictor + measure, CMP, mean)
sm$n_above <- sm$n_above * 4   # counts, over the four l_max
print(sm[order(sm$predictor, sm$measure), ], row.names=FALSE, digits=3)
cat(sprintf("\n  removals beyond the null 95th pct: %d of %d (chance ~ %.0f)\n",
    sum(CMP$n_above), 6*nrow(CMP), 0.05*6*nrow(CMP)))

## ---- impact vs module properties ---------------------------------------
PROP <- read.csv(file.path(repo,"Metrics","module_properties.csv"), stringsAsFactors=FALSE)
PROPS <- c("T_s_module","density_within","intra_density_within",
           "inter_density_within","density_across","impurity_across")
out <- list()
for (pv in names(PREDICTORS)) for (mv in MEASURES) for (lm in sort(unique(L$l_max))) {
  sub <- LX[LX$predictor==pv & LX$measure==mv & LX$l_max==lm, ]
  mm  <- merge(sub, PROP[PROP$l_max==lm, ], by.x="removed", by.y="module")
  for (p in PROPS) {
    x <- mm[[p]]; y <- abs(mm$impact)
    out[[length(out)+1]] <- data.frame(predictor=pv, measure=mv, l_max=lm, property=p,
      pearson=if (length(unique(x))>2) cor(x,y) else NA_real_,
      spearman=if (length(unique(x))>2) cor(x,y,method="spearman") else NA_real_,
      n=nrow(mm), stringsAsFactors=FALSE)
  }
}
O <- do.call(rbind, out)
write.csv(O, file.path(repo,"Metrics","loo_impact_vs_properties.csv"), row.names=FALSE)
cat("\n=== |impact| vs property: mean Spearman (n = 6 modules each)\n")
a2 <- aggregate(spearman ~ property + predictor, O, function(z) mean(z, na.rm=TRUE))
print(a2[order(a2$predictor, -abs(a2$spearman)), ], row.names=FALSE, digits=3)
cat("\nwrote loo_correlations.csv, loo_null.csv, loo_null_summary.csv, loo_impact_vs_properties.csv\n")
