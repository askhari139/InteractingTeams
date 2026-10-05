## clustering_analysis.R -------------------------------------------------
## Q3/Q4 of analysis.md: which module has the most inconsistent clustering of
## networks by phenotypic score, and which module property predicts it.
##
## The document gives no method for Q3, and "inconsistent" admits two readings.
## Both are computed, kept separate, and reported side by side:
##
##  (A) STABILITY. Cluster all 56 networks on their P_s profile. Remove module
##      X's networks, recluster the remaining 26, and compare that to the
##      full-corpus labels restricted to the same 26 (adjusted Rand index).
##      inconsistency_A = 1 - ARI. This mirrors the leave-one-out method the
##      document specifies for Q1/Q2.
##
##  (C) EXPLANATORY POWER. Adjusted Rand index between the full clustering and
##      the binary "contains X" split: how far module membership alone accounts
##      for the clustering. High = the clustering IS largely about that module.
##
##  (B) COHERENCE. On the full 56-network clustering, how scattered are the
##      networks that contain X? inconsistency_B = normalised Shannon entropy
##      of their cluster distribution (0 = all in one cluster, 1 = uniform).
##      This reads "inconsistent" as "module X's networks do not group".
##
## Features: the eight P_s columns at a given l_max, z-scored. Ward linkage on
## Euclidean distance, cut at k = 2, 3, 4.
## Null for (A): the same ARI from removing a random 26-of-56.
##
## Writes clustering_inconsistency.csv, clustering_null.csv,
## clustering_vs_properties.csv and figures.
## -----------------------------------------------------------------------
suppressMessages({library(mclust); library(cluster)})
repo <- getwd()
M   <- read.csv("Metrics/metrics_module_level.csv",  stringsAsFactors = FALSE)
NET <- read.csv("Metrics/metrics_network_level.csv", stringsAsFactors = FALSE)
PROP<- read.csv("Metrics/module_properties.csv",     stringsAsFactors = FALSE)
MODS <- c("AS","C","E","M","PS","RS")
KS <- 2:4
modsOf <- lapply(split(M$module, M$network), unique)

## feature matrix: one row per network, eight P_s columns, z-scored
features <- function(lm) {
  Ml <- M[M$l_max == lm, ]; Nl <- NET[NET$l_max == lm, ]
  agg <- function(col) as.numeric(tapply(Ml[[col]], Ml$network, mean)[Nl$network])
  X <- cbind(
    psM_infl_top30  = agg("P_s_infl_top30"),  psM_infl_top60  = agg("P_s_infl_top60"),
    psM_state_top30 = agg("P_s_state_top30"), psM_state_top60 = agg("P_s_state_top60"),
    psN_infl_top30  = Nl$P_s_net_infl_top30,  psN_infl_top60  = Nl$P_s_net_infl_top60,
    psN_state_top30 = Nl$P_s_net_state_top30, psN_state_top60 = Nl$P_s_net_state_top60)
  rownames(X) <- Nl$network
  X <- X[, apply(X, 2, function(z) sd(z) > 0), drop = FALSE]
  scale(X)
}
labels_of <- function(X, k) cutree(hclust(dist(X), method = "ward.D2"), k = k)

rows <- list(); nulls <- list(); nullsB <- list(); coh <- list()
set.seed(1); NDRAW <- 400
for (lm in sort(unique(M$l_max))) {
  X <- features(lm)
  sil <- sapply(KS, function(k) {
    l <- labels_of(X, k); mean(silhouette(l, dist(X))[, "sil_width"]) })
  for (k in KS) {
    full <- labels_of(X, k)
    ## (A) stability
    for (Xm in MODS) {
      keep <- !vapply(rownames(X), function(n) Xm %in% modsOf[[n]], logical(1))
      sub  <- labels_of(X[keep, , drop = FALSE], k)
      ari  <- adjustedRandIndex(full[keep], sub)
      rows[[length(rows)+1]] <- data.frame(l_max = lm, k = k, module = Xm,
        n_kept = sum(keep), ARI = ari, inconsistency_A = 1 - ari,
        silhouette = sil[k - min(KS) + 1], stringsAsFactors = FALSE)
    }
    ## null for (A)
    nd <- replicate(NDRAW, {
      i <- sample(nrow(X), 26)
      adjustedRandIndex(full[i], labels_of(X[i, , drop = FALSE], k)) })
    nulls[[length(nulls)+1]] <- data.frame(l_max = lm, k = k,
      draw = seq_len(NDRAW), ARI = nd, inconsistency_A = 1 - nd)
    ## null for (B): entropy of a RANDOM 30-of-56 subset's cluster spread.
    ## Needed because cluster sizes are uneven -- a module sitting entirely in
    ## the largest cluster may be no more coherent than chance.
    ndB <- replicate(NDRAW, {
      i <- sample(nrow(X), 30); tb <- table(factor(full[i], levels = 1:k))
      pp <- tb[tb > 0] / sum(tb); -sum(pp * log(pp)) / log(k) })
    nulls[[length(nulls)]]$inconsistency_B_null_med <- median(ndB)
    nullsB[[length(nullsB)+1]] <- data.frame(l_max = lm, k = k,
      draw = seq_len(NDRAW), inconsistency_B = ndB)

    ## (B) coherence of module membership in the full clustering
    for (Xm in MODS) {
      has <- vapply(rownames(X), function(n) Xm %in% modsOf[[n]], logical(1))
      tb  <- table(full[has]); p <- tb / sum(tb)
      H   <- -sum(p * log(p))
      coh[[length(coh)+1]] <- data.frame(l_max = lm, k = k, module = Xm,
        n_with = sum(has), max_frac = max(p),
        inconsistency_B = if (k > 1) H / log(k) else 0,
        membership_ARI  = adjustedRandIndex(full, has), stringsAsFactors = FALSE)
    }
  }
}
A <- do.call(rbind, rows); NU <- do.call(rbind, nulls); B <- do.call(rbind, coh)
NUB <- do.call(rbind, nullsB)
write.csv(NUB, "Metrics/clustering_null_B.csv", row.names = FALSE)
CL <- merge(A, B, by = c("l_max","k","module"))
write.csv(CL, "Metrics/clustering_inconsistency.csv", row.names = FALSE)
write.csv(NU, "Metrics/clustering_null.csv", row.names = FALSE)

cat("\n=== cluster quality (mean silhouette of the full 56-network clustering)\n")
print(reshape(unique(A[, c("l_max","k","silhouette")]), idvar = "l_max",
              timevar = "k", direction = "wide"), row.names = FALSE, digits = 3)

cat("\n=== Q3(A) stability: mean 1 - ARI when each module is removed\n")
sA <- aggregate(inconsistency_A ~ module, A, mean)
sA <- sA[order(-sA$inconsistency_A), ]
print(sA, row.names = FALSE, digits = 3)
cat(sprintf("  random-removal null: median %.3f, 5th-95th pct [%.3f, %.3f]\n",
    median(NU$inconsistency_A), quantile(NU$inconsistency_A,.05),
    quantile(NU$inconsistency_A,.95)))
above <- sum(vapply(seq_len(nrow(A)), function(i) {
  nd <- NU$inconsistency_A[NU$l_max==A$l_max[i] & NU$k==A$k[i]]
  A$inconsistency_A[i] > quantile(nd, .95) }, logical(1)))
cat(sprintf("  module removals beyond the null 95th pct: %d of %d (chance ~%.0f)\n",
    above, nrow(A), 0.05*nrow(A)))

cat("\n=== Q3(B) coherence: mean normalised entropy of cluster membership\n")
sB <- aggregate(cbind(inconsistency_B, max_frac) ~ module, B, mean)
sB <- sB[order(-sB$inconsistency_B), ]
print(sB, row.names = FALSE, digits = 3)

cat(sprintf("  random-subset null (30 of 56): median %.3f, 5th-95th pct [%.3f, %.3f]\n",
    median(NUB$inconsistency_B), quantile(NUB$inconsistency_B,.05),
    quantile(NUB$inconsistency_B,.95)))
belowB <- sum(vapply(seq_len(nrow(B)), function(i) {
  nd <- NUB$inconsistency_B[NUB$l_max==B$l_max[i] & NUB$k==B$k[i]]
  B$inconsistency_B[i] < quantile(nd, .05) }, logical(1)))
aboveB <- sum(vapply(seq_len(nrow(B)), function(i) {
  nd <- NUB$inconsistency_B[NUB$l_max==B$l_max[i] & NUB$k==B$k[i]]
  B$inconsistency_B[i] > quantile(nd, .95) }, logical(1)))
cat(sprintf("  more coherent than the null 5th pct: %d of %d; more scattered than the 95th pct: %d (chance ~%.0f each)\n",
    belowB, nrow(B), aboveB, 0.05*nrow(B)))
cat("\n=== Q3(C) explanatory power: ARI between the clustering and 'contains module X'\n")
mc <- reshape(aggregate(membership_ARI ~ module + l_max, B, mean),
              idvar = "module", timevar = "l_max", direction = "wide")
mc <- mc[order(-rowMeans(mc[, -1])), ]
print(mc, row.names = FALSE, digits = 2)

cat("\n=== cluster sizes of the full 56-network clustering\n")
for (lm in sort(unique(A$l_max))) for (k in KS) {
  X <- features(lm); full <- labels_of(X, k)
  cat(sprintf("  l_max %2d, k = %d : %s\n", lm, k,
      paste(as.integer(table(full)), collapse = " / ")))
}

cat("\n=== Q4: inconsistency vs module properties (n = 6 modules)\n")
PROPS <- c("T_s_module","density_within","intra_density_within",
           "inter_density_within","density_across","impurity_across")
out <- list()
for (lm in sort(unique(CL$l_max))) for (k in KS) {
  mm <- merge(CL[CL$l_max==lm & CL$k==k, ], PROP[PROP$l_max==lm, ], by = "module")
  for (p in PROPS) for (meas in c("inconsistency_A","inconsistency_B","membership_ARI")) {
    x <- mm[[p]]; y <- mm[[meas]]
    out[[length(out)+1]] <- data.frame(l_max=lm, k=k, property=p, measure=meas,
      spearman = if (length(unique(x))>2 && length(unique(y))>1)
                   suppressWarnings(cor(x,y,method="spearman")) else NA_real_,
      stringsAsFactors = FALSE)
  }
}
O <- do.call(rbind, out)
write.csv(O, "Metrics/clustering_vs_properties.csv", row.names = FALSE)
a2 <- aggregate(spearman ~ property + measure, O, function(z) mean(z, na.rm=TRUE))
print(a2[order(a2$measure, -abs(a2$spearman)), ], row.names = FALSE, digits = 3)
cat("\nwrote clustering_inconsistency.csv, clustering_null.csv, clustering_vs_properties.csv\n")
