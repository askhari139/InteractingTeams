## clustering_plots.R ----------------------------------------------------
## Figures for the Q3/Q4 clustering analysis.
## -----------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggrepel); library(funcsKishore)})
repo <- getwd(); outdir <- file.path(repo, "Metrics", "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

CL <- read.csv("Metrics/clustering_inconsistency.csv", stringsAsFactors = FALSE)
NU <- read.csv("Metrics/clustering_null.csv",          stringsAsFactors = FALSE)
PROP <- read.csv("Metrics/module_properties.csv",      stringsAsFactors = FALSE)

MODCOL <- c(AS="#2a78d6", C="#1baf7a", E="#eda100", M="#008300", PS="#4a3aa7", RS="#e34948")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#0b0b0b"; MUTED <- "#52514e"
thm <- function(base = 11) theme_Publication(base_size = base) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        plot.title    = element_text(hjust = 0, size = rel(1.05)),
        plot.subtitle = element_text(hjust = 0, colour = MUTED, size = rel(0.72)),
        plot.caption  = element_text(hjust = 0, colour = MUTED, size = rel(0.64)),
        axis.title    = element_text(size = rel(0.95)),
        legend.title  = element_text(size = rel(0.75)),
        legend.text   = element_text(size = rel(0.72)),
        legend.key.size = unit(0.35, "cm"),
        strip.text    = element_text(size = rel(0.8)))

## ---- 1. stability (reading A) ------------------------------------------
NU$kf <- factor(paste0("k = ", NU$k)); CL$kf <- factor(paste0("k = ", CL$k))
NU$lmf <- factor(NU$l_max); CL$lmf <- factor(CL$l_max)
p1 <- ggplot(NU, aes(lmf, inconsistency_A)) +
  geom_violin(fill = "grey88", colour = NA, width = 0.9, scale = "width") +
  geom_point(data = CL, aes(lmf, inconsistency_A), colour = BLUE, size = 1.9) +
  geom_text_repel(data = CL, aes(lmf, inconsistency_A, label = module), size = 2.4,
                  colour = INK, segment.colour = "grey70", segment.size = 0.25,
                  min.segment.length = 0.1, max.overlaps = 24, box.padding = 0.18) +
  facet_wrap(~kf, nrow = 1) +
  labs(title = "Reading A - clustering stability: how much the network clustering changes when a module is removed",
       subtitle = paste("y = 1 - adjusted Rand index between the full-corpus clustering restricted to the 26 retained networks",
                        "and the clustering recomputed on those 26. Grey violin: the same quantity for a random 26-of-56.", sep = "\n"),
       x = expression(l[max]), y = "1 - ARI",
       caption = "The null spans almost the whole range, so this reading cannot separate the modules: 3 of 72 removals beat it, chance gives ~4.") +
  thm()
ggsave(file.path(outdir, "clustering_stability.png"), p1,
       width = 11.5, height = 5.2, dpi = 200, bg = "white")

## ---- 2. coherence (reading B) ------------------------------------------
p2 <- ggplot(CL, aes(reorder(module, inconsistency_B), inconsistency_B)) +
  geom_point(aes(colour = module), size = 2.3, alpha = 0.9) +
  scale_colour_manual(values = MODCOL, guide = "none") +
  facet_wrap(~kf, nrow = 1) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Reading B - clustering coherence: how scattered a module's networks are across clusters",
       subtitle = paste("y = normalised Shannon entropy of the cluster distribution of the 30 networks containing each module,",
                        "in the full 56-network clustering. 0 = all in one cluster, 1 = spread uniformly. One point per l_max.", sep = "\n"),
       x = "module", y = "normalised entropy of cluster membership") +
  thm()
ggsave(file.path(outdir, "clustering_coherence.png"), p2,
       width = 11.5, height = 5.2, dpi = 200, bg = "white")

## ---- 3. inconsistency vs module properties -----------------------------
PROPS <- c(T_s_module = "team strength", density_within = "density within",
           intra_density_within = "intra-team density", inter_density_within = "inter-team density",
           density_across = "density across", impurity_across = "impurity across")
agg <- aggregate(cbind(inconsistency_A, inconsistency_B, membership_ARI) ~ module, CL, mean)
mm  <- merge(agg, PROP[PROP$l_max == 10, ], by = "module")
long <- do.call(rbind, lapply(names(PROPS), function(p) rbind(
  data.frame(module=mm$module, y=mm$inconsistency_A, x=mm[[p]],
             property=PROPS[[p]], reading="A: stability (1 - ARI)"),
  data.frame(module=mm$module, y=mm$inconsistency_B, x=mm[[p]],
             property=PROPS[[p]], reading="B: coherence (entropy)"),
  data.frame(module=mm$module, y=mm$membership_ARI, x=mm[[p]],
             property=PROPS[[p]], reading="C: membership explains clusters (ARI)"))))
long$property <- factor(long$property, levels = unname(PROPS))
ann <- do.call(rbind, lapply(split(long, list(long$property, long$reading), drop = TRUE),
  function(d) data.frame(property=d$property[1], reading=d$reading[1],
    lab=sprintf("rho=%+.2f", cor(d$x, d$y, method="spearman")))))
p3 <- ggplot(long, aes(x, y)) +
  geom_point(aes(colour = module), size = 2.4) +
  geom_text_repel(aes(label = module, colour = module), size = 2.2,
                  segment.colour = "grey70", segment.size = 0.22,
                  min.segment.length = 0.1, show.legend = FALSE, box.padding = 0.22) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5, size = 2.4, colour = MUTED) +
  scale_colour_manual(values = MODCOL, name = "module") +
  facet_grid(reading ~ property, scales = "free", switch = "y") +
  labs(title = "Q4 - clustering inconsistency against module properties",
       subtitle = paste("Averaged over l_max and k. One point per module; six points per panel,",
                        "so rho orients only and cannot support a claim.", sep = "\n"),
       x = "property value", y = "inconsistency") +
  thm(10) + theme(strip.text.y.left = element_text(angle = 0, hjust = 0))
ggsave(file.path(outdir, "clustering_vs_properties.png"), p3,
       width = 15.0, height = 9.4, dpi = 200, bg = "white")

## ---- 3b. reading C: does module membership explain the clustering? -----
mc <- aggregate(membership_ARI ~ module + l_max, CL, mean)
p3b <- ggplot(mc, aes(factor(l_max), membership_ARI, colour = module, group = module)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.3) +
  geom_text_repel(data = mc[mc$l_max == 10, ], aes(label = module), size = 2.6,
                  nudge_x = 0.25, segment.colour = "grey70", segment.size = 0.22,
                  show.legend = FALSE) +
  scale_colour_manual(values = MODCOL, name = "module") +
  labs(title = "Reading C - only PS membership explains the network clustering, and increasingly so with l_max",
       subtitle = paste("Adjusted Rand index between the full 56-network clustering and the binary 'contains module X' split,",
                        "averaged over k = 2, 3, 4. Zero means module membership tells you nothing about the clustering.", sep = "\n"),
       x = expression(l[max]), y = "ARI (clustering vs contains-module)") +
  thm()
ggsave(file.path(outdir, "clustering_membership_ARI.png"), p3b,
       width = 9.0, height = 5.2, dpi = 200, bg = "white")

## ---- 4. what the clusters actually look like ---------------------------
## cluster composition per module, at l_max = 10, k = 3
suppressMessages({library(mclust)})
M   <- read.csv("Metrics/metrics_module_level.csv",  stringsAsFactors = FALSE)
NET <- read.csv("Metrics/metrics_network_level.csv", stringsAsFactors = FALSE)
lm <- 10; k <- 3
Ml <- M[M$l_max == lm, ]; Nl <- NET[NET$l_max == lm, ]
agg2 <- function(col) as.numeric(tapply(Ml[[col]], Ml$network, mean)[Nl$network])
X <- cbind(agg2("P_s_infl_top30"), agg2("P_s_infl_top60"),
           agg2("P_s_state_top30"), agg2("P_s_state_top60"),
           Nl$P_s_net_infl_top30, Nl$P_s_net_infl_top60,
           Nl$P_s_net_state_top30, Nl$P_s_net_state_top60)
rownames(X) <- Nl$network; X <- scale(X)
lab <- cutree(hclust(dist(X), method = "ward.D2"), k = k)
modsOf <- lapply(split(M$module, M$network), unique)
comp <- do.call(rbind, lapply(names(MODCOL), function(m) {
  has <- vapply(rownames(X), function(n) m %in% modsOf[[n]], logical(1))
  tb <- table(factor(lab[has], levels = 1:k))
  data.frame(module = m, cluster = factor(1:k), frac = as.numeric(tb/sum(tb)),
             n = as.numeric(tb)) }))
p4 <- ggplot(comp, aes(cluster, frac, fill = module)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = ifelse(n > 0, n, "")), vjust = -0.4, size = 2.6, colour = INK) +
  scale_fill_manual(values = MODCOL, guide = "none") +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~module, nrow = 1) +
  labs(title = sprintf("How each module's 30 networks distribute over %d clusters (l_max = %d)", k, lm),
       subtitle = paste("Clusters from Ward linkage on the eight z-scored P_s columns, all 56 networks.",
                        "Labels are network counts. PS concentrates in one cluster; M and C spread.", sep = "\n"),
       x = "cluster", y = "fraction of the module's networks") +
  thm()
ggsave(file.path(outdir, "clustering_composition.png"), p4,
       width = 12.0, height = 4.6, dpi = 200, bg = "white")

cat("wrote 4 clustering figures\n")
