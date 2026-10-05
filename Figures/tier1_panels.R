## tier1_panels.R ---------------------------------------------------------
## Tier-1 panels of FIGURE_PLAN.md: every panel whose inputs are already on
## disk and validated. Writes Figures/panels/*.png plus two composites.
##
##   fig1b   module property heatmap (impurity_within dropped, see plan)
##   fig1h   module T_s vs module P_s, dominant and top60
##   fig1i   network T_s vs PC1 variance
##   fig2    6x6 module-pair heatmaps, 4 tile variables            [composite]
##   fig3d   directed inter-module asymmetry by module count
##   fig4a   network T_s vs network P_s            (the r = 0.94 anchor)
##   fig4b   inter-module T_s sum vs the same P_s  (the sign contrast)
##   fig4d   leave-one-out impact by module, against the resampling null
##   fig4e   leave-one-out violins
##   fig4f   clustering reading C: membership ARI by l_max
##
## Conventions (FIGURE_PLAN.md): l_max = 10, influence-weighted P_s, top60.
## Run from the repo root:  Rscript Figures/tier1_panels.R
## -------------------------------------------------------------------------
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggrepel)
  library(patchwork); library(funcsKishore)
})

outdir <- file.path("Figures", "panels")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
LMAX <- 10

## --- palette -------------------------------------------------------------
## Modules: the 6-of-8 subset of the reference categorical palette already
## validated for this repo (Metrics/loo_plots.R). Worst-pair CVD dE 6.9 sits
## in the band that requires secondary encoding, so every module mark is also
## direct-labelled.
MODCOL <- c(AS="#2a78d6", C="#1baf7a", E="#eda100", M="#008300",
            PS="#4a3aa7", RS="#e34948")
## Module count is ordinal -> one hue, light->dark, no step lighter than 250.
NMODCOL <- c(`2`="#86b6ef", `3`="#5598e7", `4`="#2a78d6", `5`="#184f95")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; RED <- "#e34948"
GREY <- "#f0efec"; INK <- "#0b0b0b"; MUTED <- "#52514e"
SEQ  <- c("#cde2fb","#9ec5f4","#6da7ec","#3987e5","#256abf","#184f95","#0d366b")

## --- provenance ----------------------------------------------------------
## Ising is the simulation formalism for this paper: every phenotypic score is
## computed from the Boolean/Ising steady states in <net>_finFlagFreq.csv.
## RACIPE contributes the PCA spectrum ONLY (<Tag>_PCA_variance.csv, 58k-265k
## continuous samples per network). Topology and team-strength metrics come
## from the influence matrix and involve no simulation at all.
SRC_ISING  <- "Phenotypic scores: Ising steady states (finFlagFreq, basin-weighted)."
SRC_RACIPE <- "PCA spectrum: RACIPE continuous sampling (58k-265k samples per network)."
SRC_TOPO   <- "Structural metrics from the influence matrix; no simulation."
## wrap captions: an unwrapped one-liner runs off the canvas
cap <- function(..., width = 110) paste(strwrap(paste(c(...), collapse = " "), width = width), collapse = "\n")

thm <- function(base = 11) theme_Publication(base_size = base) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        plot.title    = element_text(hjust = 0, size = rel(1.0), face = "bold"),
        plot.subtitle = element_text(hjust = 0, colour = MUTED, size = rel(0.72)),
        plot.caption  = element_text(hjust = 0, colour = MUTED, size = rel(0.62)),
        axis.title    = element_text(size = rel(0.92)),
        legend.title  = element_text(size = rel(0.75)),
        legend.text   = element_text(size = rel(0.72)),
        legend.key.size = unit(0.35, "cm"),
        strip.text    = element_text(size = rel(0.82)))
thm_tile <- function(base = 11) thm(base) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank())

sv <- function(p, f, w, h) {
  ggsave(file.path(outdir, f), p, width = w, height = h, dpi = 400, bg = "white")
  cat(sprintf("  %-34s %4.1f x %4.1f in\n", f, w, h))
}
rlab <- function(x, y) sprintf("r = %+.3f", cor(x, y))

## --- data ----------------------------------------------------------------
N  <- read.csv("Metrics/metrics_network_level.csv", stringsAsFactors = FALSE)
M  <- read.csv("Metrics/metrics_module_level.csv",  stringsAsFactors = FALSE)
P  <- read.csv("Metrics/metrics_pair_level.csv",    stringsAsFactors = FALSE)
PR <- read.csv("Metrics/module_properties.csv",     stringsAsFactors = FALSE)
L  <- read.csv("Metrics/loo_correlations.csv",      stringsAsFactors = FALSE)
NU <- read.csv("Metrics/loo_null.csv",              stringsAsFactors = FALSE)
CI <- read.csv("Metrics/clustering_inconsistency.csv", stringsAsFactors = FALSE)

n10 <- N[N$l_max == LMAX, ]; m10 <- M[M$l_max == LMAX, ]
p10 <- P[P$l_max == LMAX, ]; pr10 <- PR[PR$l_max == LMAX, ]
n10$nmod <- factor(n10$n_modules)

## module order: weakest to strongest team, at l_max 10
MODORD <- pr10$module[order(pr10$T_s_module)]
cat("module order by T_s (l_max 10): ", paste(MODORD, collapse = " < "), "\n\n")

cat("writing panels to ", outdir, "\n", sep = "")

## =========================================================================
## Fig 1b -- module property heatmap
## impurity_within is dropped: it is 0 for 5 of 6 modules (PS alone, 0.143),
## so as a heatmap row it is a blank strip. It stays in the S1 table.
## =========================================================================
PROPS <- c(T_s_module          = "team strength",
           density_within      = "density (within)",
           intra_density_within= "intra-team density",
           inter_density_within= "inter-team density",
           density_across      = "density (across)",
           impurity_across     = "impurity (across)")

hm <- pr10 %>%
  select(module, all_of(names(PROPS))) %>%
  pivot_longer(-module, names_to = "prop", values_to = "raw") %>%
  group_by(prop) %>% mutate(z = as.numeric(scale(raw))) %>% ungroup() %>%
  mutate(module = factor(module, levels = MODORD),
         prop   = factor(PROPS[prop], levels = PROPS),
         lab    = ifelse(raw < 0.02, sprintf("%.3f", raw), sprintf("%.2f", raw)))

p1b <- ggplot(hm, aes(module, prop, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = lab, colour = abs(z) > 1.1), size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                       name = "z-score\n(within row)") +
  scale_y_discrete(limits = rev(levels(hm$prop))) +
  labs(title = "Module properties",
       subtitle = "six modules, ordered by team strength; tiles scaled within each property, cells give the raw value",
       x = NULL, y = NULL,
       caption = cap("l_max = 10.", SRC_TOPO, "impurity (within) omitted: 0 for every module but PS (0.143).")) +
  thm_tile() + theme(legend.position = "right")
sv(p1b, "fig1b_module_properties.png", 7.2, 3.4)

## =========================================================================
## Fig 1h -- module team strength vs module phenotypic score
## =========================================================================
d1h <- m10 %>%
  select(network, module, T_module_intra, P_s_infl_dominant, P_s_infl_top60) %>%
  pivot_longer(starts_with("P_s_"), names_to = "measure", values_to = "P_s") %>%
  mutate(measure = recode(measure,
                          P_s_infl_dominant = "dominant state",
                          P_s_infl_top60    = "top 60% of the state space"))
d1h$measure <- factor(d1h$measure,
                      levels = c("dominant state", "top 60% of the state space"))
## r within each module, to show how much of the trend is between-module
within <- d1h %>% group_by(measure, module) %>%
  summarise(rw = cor(T_module_intra, P_s), .groups = "drop") %>%
  group_by(measure) %>% summarise(mn = mean(rw), .groups = "drop")
ann1h <- d1h %>% group_by(measure) %>%
  summarise(lab = rlab(T_module_intra, P_s), .groups = "drop") %>%
  left_join(within, by = "measure") %>%
  mutate(lab = sprintf("overall %s\nmean within-module r = %+.2f", lab, mn))
## module means: the six points the overall r is really carried by
mu1h <- d1h %>% group_by(measure, module) %>%
  summarise(x = mean(T_module_intra), y = mean(P_s), .groups = "drop")

p1h <- ggplot(d1h, aes(T_module_intra, P_s)) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
              linewidth = 0.8, alpha = 0.3) +
  geom_point(aes(colour = module), size = 1.4, alpha = 0.5, stroke = 0) +
  geom_point(data = mu1h, aes(x, y, fill = module), shape = 21, size = 3.2,
             colour = "white", stroke = 0.8) +
  geom_text_repel(data = mu1h, aes(x, y, label = module, colour = module),
                  size = 2.9, fontface = "bold", seed = 1,
                  min.segment.length = 0.3, segment.size = 0.25,
                  box.padding = 0.35, show.legend = FALSE) +
  geom_text(data = ann1h, aes(x = 0.04, y = 1.62, label = lab), hjust = 0,
            vjust = 1, size = 2.8, colour = MUTED, lineheight = 1.05,
            inherit.aes = FALSE) +
  facet_wrap(~measure) +
  scale_colour_manual(values = MODCOL, guide = "none") +
  scale_fill_manual(values = MODCOL, guide = "none") +
  labs(title = "Stronger module teams give a stronger phenotypic score",
       subtitle = "180 module instances across the 56 networks; large points are the six module means",
       x = expression("module team strength  "*T[s]^{mod}),
       y = expression("module "*P[s]*"  (influence-weighted)"),
       caption = paste("l_max = 10. Read the overall r with care: module team strength is",
                       "clustered by module, so the trend is carried by differences\nbetween the six modules,",
                       "not by variation within them. Every module is direct-labelled.\n", SRC_ISING)) +
  thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(p1h, "fig1h_Ts_vs_Ps.png", 8.0, 4.4)

## =========================================================================
## Fig 1i -- team strength vs PC1 variance
## DEPARTURE from FIGURE_PLAN.md: plotted at the NETWORK level, not module.
## PC1_var is a network property, so joining it onto the 180 module rows
## would repeat each y three times over and inflate n. 56 clean points.
## =========================================================================
p1i <- ggplot(n10, aes(T_network, PC1_var)) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
              linewidth = 0.8, alpha = 0.3) +
  geom_point(aes(fill = nmod), shape = 21, size = 2.6, colour = "white", stroke = 0.7) +
  annotate("text", x = min(n10$T_network), y = max(n10$PC1_var),
           label = rlab(n10$T_network, n10$PC1_var), hjust = 0, size = 3.2) +
  scale_fill_manual(values = NMODCOL, name = "modules") +
  labs(title = "Stronger teams give a more concentrated landscape",
       subtitle = "variance explained by PC1, 56 networks",
       x = expression("network team strength  "*T[s]), y = "PC1 variance explained",
       caption = cap("l_max = 10.", SRC_RACIPE, "Team strength is structural (influence matrix).")) +
  thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(p1i, "fig1i_Ts_vs_PC1.png", 5.6, 4.0)

## =========================================================================
## Fig 2 -- 6x6 module-pair heatmaps (the 15 two-module networks)
## =========================================================================
two <- n10[n10$n_modules == 2, ]
two$mods <- lapply(strsplit(sub("^Comb_", "", two$network), "_"), sort)
two$A <- vapply(two$mods, `[`, "", 1); two$B <- vapply(two$mods, `[`, "", 2)

TILES <- c(P_s_net_infl_dominant = "P_s, dominant state  (Ising)",
           P_s_net_infl_top60    = "P_s, top 60%  (Ising)",
           PC1_var               = "PC1 variance  (RACIPE)",
           PCA_entropy_norm      = "PCA entropy, normalised  (RACIPE)")

mk2 <- function(col) {
  d <- two[, c("A", "B", col)]; names(d)[3] <- "v"
  ## the matrix is symmetric -- one network per unordered pair
  d <- rbind(d, setNames(data.frame(d$B, d$A, d$v), names(d)))
  grid <- expand.grid(A = MODORD, B = MODORD, stringsAsFactors = FALSE)
  d <- merge(grid, d, all.x = TRUE)
  d$A <- factor(d$A, levels = MODORD); d$B <- factor(d$B, levels = rev(MODORD))
  rng <- range(d$v, na.rm = TRUE)
  ## white ink only on the dark end of the ramp, by scaled position
  d$sc  <- (d$v - rng[1]) / diff(rng)
  d$ink <- !is.na(d$sc) & d$sc > 0.55
  ggplot(d, aes(A, B, fill = v)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = ifelse(is.na(v), "", sprintf("%.2f", v)), colour = ink),
              size = 2.7, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), na.value = INK) +
    scale_fill_gradientn(colours = SEQ, na.value = "grey96", name = NULL,
                         breaks = scales::breaks_pretty(3),
                         guide = guide_colourbar(barwidth  = unit(0.30, "cm"),
                                                 barheight = unit(2.4, "cm"),
                                                 ticks.colour = "white")) +
    labs(title = TILES[[col]], x = NULL, y = NULL) +
    thm_tile(10) +
    theme(legend.position = "right", legend.direction = "vertical",
          plot.title = element_text(size = rel(0.95), face = "bold"))
}
p2 <- wrap_plots(lapply(names(TILES), mk2), ncol = 2) +
  plot_annotation(
    title = "Which two modules are combined sets the landscape",
    subtitle = paste("the 15 two-module networks; both axes ordered by module team strength",
                     paste(MODORD, collapse = " < ")),
    caption = paste("l_max = 10, influence-weighted P_s. Diagonal blank: single-module networks are not in this corpus.",
                    "\nb, c:", SRC_ISING, "  d, e:", SRC_RACIPE),
    tag_levels = list(c("b", "c", "d", "e")),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(colour = MUTED, size = 9),
                  plot.caption = element_text(colour = MUTED, size = 8, hjust = 0)))
sv(p2, "fig2_pair_heatmaps.png", 9.0, 7.2)

## =========================================================================
## Fig 3d -- directed asymmetry of inter-module interaction
##   asym = |T_A->B - T_B->A| / (T_A->B + T_B->A)   in [0, 1]
## Each unordered pair appears twice in the table; keep one ordering.
## =========================================================================
pu    <- p10[p10$module_from < p10$module_to, ]
npair <- nrow(pu); ndrop <- sum(pu$T_from_to + pu$T_to_from == 0)
d3 <- pu %>%
  mutate(tot = T_from_to + T_to_from,
         asym = ifelse(tot > 0, abs(T_from_to - T_to_from) / tot, NA_real_),
         nmod = factor(n_modules)) %>%
  filter(!is.na(asym))
med3 <- d3 %>% group_by(nmod) %>%
  summarise(m = median(asym), mu = mean(asym), n = n(), .groups = "drop")

p3d <- ggplot(d3, aes(nmod, asym)) +
  geom_violin(aes(fill = nmod), colour = NA, alpha = 0.55, width = 0.85, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               colour = INK, linewidth = 0.4) +
  geom_text(data = med3, aes(nmod, 1.06, label = sprintf("mean %.2f\nn = %d", mu, n)),
            size = 2.7, colour = MUTED, lineheight = 0.95) +
  scale_fill_manual(values = NMODCOL, guide = "none") +
  scale_y_continuous(limits = c(0, 1.14), breaks = seq(0, 1, 0.25)) +
  labs(title = "Inter-module coupling is mostly one-way, less so in bigger networks",
       subtitle = expression("directed asymmetry  "*group("|",T[A%->%B]-T[B%->%A],"|")/(T[A%->%B]+T[B%->%A])*"  per ordered module pair"),
       x = "modules in the network", y = "directed asymmetry",
       caption = cap(sprintf("l_max = 10. 1 = strictly one-directional, 0 = perfectly reciprocal. %d of %d module pairs have no coupling in either direction and are dropped.",
                             ndrop, npair), "\n", SRC_TOPO)) +
  thm()
sv(p3d, "fig3d_asymmetry.png", 5.8, 4.2)

## =========================================================================
## Fig 4a / 4b -- the two anchor scatters
## =========================================================================
anchor <- function(xcol, xlab, title, sub) {
  x <- n10[[xcol]]; y <- n10$P_s_net_infl_top60
  ggplot(n10, aes(.data[[xcol]], P_s_net_infl_top60)) +
    geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
                linewidth = 0.8, alpha = 0.3) +
    geom_point(aes(fill = nmod), shape = 21, size = 2.6, colour = "white", stroke = 0.7) +
    annotate("label", x = min(x), y = max(y), label = rlab(x, y), hjust = 0, vjust = 1,
             size = 3.4, label.size = 0, fill = "white", alpha = 0.8) +
    scale_fill_manual(values = NMODCOL, name = "modules") +
    labs(title = title, subtitle = sub, x = xlab,
         y = expression("network "*P[s]*"  (influence-weighted, top 60%)"),
         caption = SRC_ISING) +
    thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
}
p4a <- anchor("T_network", expression("network team strength  "*T[s]),
              "Team strength predicts the phenotypic score",
              "56 networks, l_max = 10")
sv(p4a, "fig4a_Tnetwork_anchor.png", 5.6, 4.2)

## 4b carries an inset, because at l_max 10 this correlation is at its weakest
## (-0.21) and the l_max trend is the substantive part of the result.
trend <- do.call(rbind, lapply(c(1, 5, 7, 10), function(l) {
  d <- N[N$l_max == l, ]
  data.frame(l_max = l, r = cor(d$T_inter_sum, d$P_s_net_infl_top60))
}))
p4b <- anchor("T_inter_sum", expression("inter-module coupling  "*sum(T[A%->%B])),
              "Inter-module coupling runs the other way",
              "the sign is opposite to 4a; 56 networks, l_max = 10")
sv(p4b, "fig4b_Tinter_anchor.png", 5.6, 4.2)

## 4c: the l_max dependence, which is the substantive part -- at l_max 10 this
## correlation is at its weakest. A separate panel, not an inset: as an inset
## it collapses against the scatter it sits on.
p4c <- ggplot(trend, aes(factor(l_max), r, group = 1)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.4) +
  geom_line(colour = RED, linewidth = 1.0) +
  geom_point(colour = RED, size = 2.6) +
  geom_text(aes(label = sprintf("%+.2f", r)), vjust = 1.9, size = 3.0, colour = MUTED) +
  scale_y_continuous(limits = c(-0.68, 0.05)) +
  labs(title = "The coupling effect fades with path length",
       subtitle = "same correlation as 4b, recomputed at each path-length cutoff",
       x = expression(l[max]), y = expression("r  ( "*sum(T[A%->%B])*" vs network "*P[s]*" )"),
       caption = cap("l_max = 10, used everywhere else, is where this relationship is weakest.", SRC_ISING)) +
  thm() + theme(panel.grid.major.x = element_blank())
sv(p4c, "fig4c_lmax_trend.png", 4.6, 4.2)

## =========================================================================
## Fig 4d -- LOO impact by module, against the resampling null
## Fig 4e -- LOO violins
##   Both on the pairing that carries real signal: inter-module T_s sum
##   against the network-level influence-weighted P_s (loo_results.md).
## =========================================================================
PV <- "T_inter_sum"; FM <- "psN_infl"
lx <- L[L$predictor == PV & L$family == FM & L$removed != "All", ]
la <- L[L$predictor == PV & L$family == FM & L$removed == "All", ]
nd <- NU[NU$predictor == PV & NU$family == FM, ]
ML <- c(psN_infl_top30 = "top 30%", psN_infl_top60 = "top 60%")
lx$panel <- ML[lx$measure]; la$panel <- ML[la$measure]; nd$panel <- ML[nd$measure]

## Null band: the pooled 95th percentile of |null deviation|, the same
## criterion loo_results.md uses, so figure and text agree.
P95 <- unname(quantile(abs(nd$r - nd$r_all), 0.95))
band <- data.frame(panel = unique(lx$panel), lo = -P95, hi = P95)
nbeyond <- sum(abs(lx$impact) > P95)
modbeyond <- sort(unique(lx$removed[abs(lx$impact) > P95]))
cat(sprintf("  [4d] null p95 = %.3f; %d of %d removals clear it: %s\n",
            P95, nbeyond, nrow(lx), paste(modbeyond, collapse = ", ")))

p4d <- ggplot(lx, aes(reorder(removed, -abs(impact)), impact)) +
  geom_rect(data = band, inherit.aes = FALSE,
            aes(xmin = -Inf, xmax = Inf, ymin = lo, ymax = hi),
            fill = "grey88", alpha = 0.7) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  geom_point(aes(fill = factor(l_max)), shape = 21, size = 2.5,
             colour = "white", stroke = 0.6,
             position = position_dodge(width = 0.55)) +
  facet_wrap(~panel, nrow = 1) +
  scale_fill_manual(values = unname(NMODCOL), name = expression(l[max])) +
  labs(title = "Removing M moves this correlation most; AS is the only other module to clear the null",
       subtitle = sprintf("change in r when one module's 30 networks are dropped; grey band is the resampling null (|\u0394r| < %.2f, 95th pct of 2000 random 26-of-56 subsets)", P95),
       x = "module removed", y = expression(Delta*r),
       caption = cap(sprintf("Predictor: inter-module T_s sum (influence matrix). Outcome: network-level influence-weighted P_s. %d of %d removals clear the null.", nbeyond, nrow(lx)), "\n", SRC_ISING)) +
  thm() + theme(panel.grid.major.x = element_blank())
sv(p4d, "fig4d_loo_impact.png", 7.0, 4.0)

lx$kind <- "one module removed"; la$kind <- "all 56 networks"
p4e <- ggplot(nd, aes(factor(l_max), r)) +
  geom_violin(fill = "grey88", colour = NA, scale = "width", width = 0.9) +
  geom_point(data = lx, aes(factor(l_max), r), colour = BLUE, size = 1.8) +
  geom_text_repel(data = lx, aes(factor(l_max), r, label = removed),
                  size = 2.3, colour = MUTED, min.segment.length = 0.2,
                  segment.size = 0.25, max.overlaps = 20, box.padding = 0.18) +
  geom_point(data = la, aes(factor(l_max), r), colour = ORANGE, shape = 18, size = 3.4) +
  facet_wrap(~panel, nrow = 1) +
  labs(title = "The six module removals sit inside the resampling spread",
       subtitle = "grey violin: r over 2000 random 26-of-56 subsets  |  blue: the six leave-one-module-out values  |  orange diamond: all 56 networks",
       x = expression(l[max]), y = "Pearson r",
       caption = cap("Predictor: inter-module T_s sum (influence matrix). Outcome: network-level influence-weighted P_s.", "\n", SRC_ISING)) +
  thm() + theme(panel.grid.major.x = element_blank())
sv(p4e, "fig4e_loo_violins.png", 7.6, 4.2)

## =========================================================================
## Fig 4f -- clustering reading C: does membership explain the clustering?
## =========================================================================
d4f <- CI %>% group_by(module, l_max) %>%
  summarise(ARI = mean(membership_ARI), .groups = "drop")
lastpt <- d4f[d4f$l_max == max(d4f$l_max), ]

p4f <- ggplot(d4f, aes(l_max, ARI, colour = module, group = module)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_text_repel(data = lastpt, aes(label = module), hjust = 0, nudge_x = 0.5,
                  direction = "y", size = 3.1, segment.size = 0.25,
                  min.segment.length = 0.3, show.legend = FALSE) +
  scale_colour_manual(values = MODCOL, name = "module") +
  scale_x_continuous(breaks = c(1, 5, 7, 10), limits = c(1, 12.2)) +
  labs(title = "PS, and only PS, explains how the networks cluster",
       subtitle = "ARI between the phenotypic-score clustering and the binary \"contains this module\" split, averaged over k = 2, 3, 4",
       x = expression(l[max]), y = "adjusted Rand index",
       caption = cap("ARI 0 = membership tells you nothing about the clustering. Clustering features are the eight P_s columns.", "\n", SRC_ISING)) +
  thm() + theme(legend.position = "none",
                panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(p4f, "fig4f_clustering_ARI.png", 5.8, 4.2)

## =========================================================================
## Composite: Figure 4
## =========================================================================
fig4 <- (p4a | p4b | p4c) / (p4d | p4f) +
  plot_annotation(tag_levels = "a",
    title = "Figure 4  Do the interaction metrics predict the landscape?",
    theme = theme(plot.title = element_text(face = "bold", size = 14)))
sv(fig4, "FIG4_composite.png", 16.0, 9.0)

## Composite: Figure 1 (the panels that exist so far)
fig1 <- (p1b) / (p1h | p1i) +
  plot_layout(heights = c(0.85, 1.15)) +
  plot_annotation(tag_levels = list(c("b", "h", "i")),
    title = "Figure 1  Two-team modules and their landscapes (Tier-1 panels only)",
    caption = "Panels a, c-g of FIGURE_PLAN.md are not buildable from data now on disk.",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.caption = element_text(colour = MUTED, size = 8, hjust = 0)))
sv(fig1, "FIG1_partial_composite.png", 13.0, 9.5)

cat("\ndone.\n")
