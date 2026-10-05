## loo_plots.R ------------------------------------------------------------
## Figures for the leave-one-module-out analysis, per predictor and per
## outcome family. Uses funcsKishore::theme_Publication.
##
## Per predictor (network T_s, inter-module T_s sum):
##   violin_<predictor>_<family>.png   distribution of correlations: the
##       random-subset null as a violin, the six leave-one-module-out values,
##       and the original all-56-network correlation
##   impact_<predictor>.png            change in r by module, coloured by the
##       exact correlation the removal produced
##   properties_<predictor>.png        change in r against module properties,
##       coloured by module
##
## Families are plotted separately because their correlation scales differ by
## an order of magnitude -- pooling them would hide that.
## -------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(ggrepel); library(funcsKishore)})
repo <- getwd(); outdir <- file.path(repo, "Metrics", "figures")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
unlink(Sys.glob(file.path(outdir, "*.png")))

L     <- read.csv("Metrics/loo_correlations.csv", stringsAsFactors = FALSE)
NULLD <- read.csv("Metrics/loo_null.csv",         stringsAsFactors = FALSE)
PROP  <- read.csv("Metrics/module_properties.csv",stringsAsFactors = FALSE)

## --- palette -------------------------------------------------------------
## Module hues: the 6-of-8 subset of the reference categorical palette that
## clears the all-pairs colourblind (6.9) and normal-vision (15.6) floors.
## CVD in the 6-8 band is legal only with secondary encoding, so every module
## point also carries a direct text label.
MODCOL <- c(AS="#2a78d6", C="#1baf7a", E="#eda100", M="#008300",
            PS="#4a3aa7", RS="#e34948")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; RED <- "#e34948"
GREY <- "#f0efec";  INK <- "#0b0b0b";   MUTED <- "#52514e"

PRED <- c(T_network = "network T_s", T_inter_sum = "inter-module T_s sum")
FAMLAB <- c(psM_infl   = "module-level P_s, influence-weighted",
            psM_state  = "module-level P_s, state-based",
            psN_infl   = "network-level P_s, influence-weighted",
            psN_state  = "network-level P_s, state-based",
            pc         = "PCA metrics")
## short labels, used inside a family where they are unique
MLAB <- c(psM_infl_top30="top30", psM_infl_top60="top60",
          psM_state_top30="top30", psM_state_top60="top60",
          psN_infl_top30="top30", psN_infl_top60="top60",
          psN_state_top30="top30", psN_state_top60="top60",
          pc_PC1_var="PC1 variance", pc_NumPC_90="PCs to 90%",
          pc_PCA_entropy="PCA entropy", pc_PCA_entropy_norm="PCA entropy (norm)")
## unique labels, used when faceting across families
PANEL <- c(psM_infl_top30="P_s module, infl, top30", psM_infl_top60="P_s module, infl, top60",
           psM_state_top30="P_s module, state, top30", psM_state_top60="P_s module, state, top60",
           psN_infl_top30="P_s network, infl, top30", psN_infl_top60="P_s network, infl, top60",
           psN_state_top30="P_s network, state, top30", psN_state_top60="P_s network, state, top60",
           pc_PC1_var="PC1 variance", pc_NumPC_90="PCs to 90%",
           pc_PCA_entropy="PCA entropy", pc_PCA_entropy_norm="PCA entropy (norm)")

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

LX <- L[L$removed != "All", ]; LA <- L[L$removed == "All", ]

## ---- 1. violins, per predictor x family --------------------------------
for (pv in names(PRED)) for (fm in names(FAMLAB)) {
  nd <- NULLD[NULLD$predictor == pv & NULLD$family == fm, ]
  lx <- LX[LX$predictor == pv & LX$family == fm, ]
  la <- LA[LA$predictor == pv & LA$family == fm, ]
  if (!nrow(nd)) next
  for (D in list(nd, lx, la)) NULL
  nd$panel <- MLAB[nd$measure]; lx$panel <- MLAB[lx$measure]; la$panel <- MLAB[la$measure]
  nd$lmf <- factor(nd$l_max); lx$lmf <- factor(lx$l_max); la$lmf <- factor(la$l_max)
  la$kind <- "all 56 networks (original)"; lx$kind <- "one module removed"

  p <- ggplot(nd, aes(lmf, r)) +
    geom_violin(fill = "grey88", colour = NA, width = 0.9, scale = "width") +
    geom_point(data = lx, aes(lmf, r, colour = kind, shape = kind), size = 1.9) +
    geom_point(data = la, aes(lmf, r, colour = kind, shape = kind), size = 3.1) +
    geom_text_repel(data = lx, aes(lmf, r, label = removed), size = 2.4,
                    colour = INK, segment.colour = "grey70", segment.size = 0.25,
                    min.segment.length = 0.1, max.overlaps = 24, box.padding = 0.18) +
    scale_colour_manual(values = c("one module removed" = BLUE,
                                   "all 56 networks (original)" = ORANGE), name = NULL) +
    scale_shape_manual(values = c("one module removed" = 16,
                                  "all 56 networks (original)" = 18), name = NULL) +
    facet_wrap(~panel, scales = "free_y", nrow = 1) +
    labs(title = sprintf("%s vs %s: correlation when one module is removed", FAMLAB[[fm]], PRED[[pv]]),
         subtitle = paste("Grey violin: r over 2000 random 26-of-56 network subsets, the same number a module removal leaves.",
                          "Blue: r with one module's networks removed, labelled by module. Orange diamond: the original r over all 56 networks.", sep = "\n"),
         x = expression(l[max]), y = "Pearson r",
         caption = "A module's effect is only interpretable against the grey violin: r on 26 networks is noisy regardless of which 26.") +
    thm()
  ggsave(file.path(outdir, sprintf("violin_%s_%s.png", pv, fm)), p,
         width = 3.1 * length(unique(nd$panel)) + 2.4, height = 5.0, dpi = 200, bg = "white")
}

## ---- 2. impact by module, coloured by the exact correlation -------------
## Split into module-level P_s, network-level P_s and PCA metrics: their
## correlation scales differ enough that a shared figure hides the contrast.
GROUP <- c(psM_infl = "psM", psM_state = "psM",
           psN_infl = "psN", psN_state = "psN", pc = "pc")
GRPLAB <- c(psM = "module-level P_s", psN = "network-level P_s", pc = "PCA metrics")
for (pv in names(PRED)) for (gr in names(GRPLAB)) {
  lx <- LX[LX$predictor == pv & GROUP[LX$family] == gr, ]
  if (!nrow(lx)) next
  lx$panel <- factor(PANEL[lx$measure], levels = unname(PANEL))
  nsel <- NULLD$predictor == pv & GROUP[NULLD$family] == gr
  band <- aggregate(list(p95 = abs(NULLD$r_all - NULLD$r)[nsel]),
                    by = list(measure = NULLD$measure[nsel]),
                    FUN = function(z) quantile(z, .95))
  band$panel <- factor(PANEL[band$measure], levels = unname(PANEL))

  p <- ggplot(lx, aes(removed, impact)) +
    geom_rect(data = band, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = -p95, ymax = p95), fill = "grey91") +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    ## shape 21 so near-zero r (pale fill) keeps an outline against the band
    geom_point(aes(fill = r), shape = 21, size = 2.5, colour = "grey35", stroke = 0.3) +
    scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                         limits = c(-1, 1), name = "r after removal") +
    facet_wrap(~panel, scales = "free_y", nrow = 1) +
    labs(title = sprintf("%s: change in the correlation with %s when each module is removed",
                         GRPLAB[[gr]], PRED[[pv]]),
         subtitle = paste("One point per l_max (1, 5, 7, 10). Positive means removing that module lowers r.",
                          "Point colour is the correlation the removal actually produced, on a fixed -1 to +1 scale.",
                          "Grey band: central 90% of the random-subset null.", sep = "\n"),
         x = "module removed", y = expression(r[all] - r[without~module])) +
    thm(10)
  ggsave(file.path(outdir, sprintf("impact_%s_%s.png", pv, gr)), p,
         width = 3.2 * length(unique(lx$panel)) + 2.6, height = 5.2, dpi = 200, bg = "white")
}

## ---- 3. impact vs module properties, coloured by module ----------------
LMF <- 10
PROPS <- c(T_s_module = "team strength", density_within = "density within",
           intra_density_within = "intra-team density", inter_density_within = "inter-team density",
           density_across = "density across", impurity_across = "impurity across")
## Split by group as well, so each figure is at most two family rows wide.
for (pv in names(PRED)) for (gr in names(GRPLAB)) {
  sub <- LX[LX$predictor == pv & LX$l_max == LMF & GROUP[LX$family] == gr, ]
  if (!nrow(sub)) next
  imp <- aggregate(list(y = abs(sub$impact)),
                   by = list(module = sub$removed, family = sub$family), FUN = mean)
  mm  <- merge(imp, PROP[PROP$l_max == LMF, ], by = "module")
  long <- do.call(rbind, lapply(names(PROPS), function(p) data.frame(
    module = mm$module, family = mm$family, y = mm$y, x = mm[[p]],
    property = PROPS[[p]], stringsAsFactors = FALSE)))
  long$property <- factor(long$property, levels = unname(PROPS))
  long$famlab   <- factor(FAMLAB[long$family], levels = unname(FAMLAB))
  ann <- do.call(rbind, lapply(split(long, list(long$property, long$famlab), drop = TRUE),
    function(d) data.frame(property = d$property[1], famlab = d$famlab[1],
      lab = sprintf("rho=%+.2f", cor(d$x, d$y, method = "spearman")))))

  p <- ggplot(long, aes(x, y)) +
    geom_point(aes(colour = module), size = 2.4) +
    geom_text_repel(aes(label = module, colour = module), size = 2.2,
                    segment.colour = "grey70", segment.size = 0.22,
                    min.segment.length = 0.1, show.legend = FALSE, box.padding = 0.22) +
    geom_text(data = ann, aes(x = Inf, y = Inf, label = lab), inherit.aes = FALSE,
              hjust = 1.1, vjust = 1.5, size = 2.4, colour = MUTED) +
    scale_colour_manual(values = MODCOL, name = "module") +
    facet_grid(famlab ~ property, scales = "free", switch = "y") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.05, 0.25))) +
    labs(title = sprintf("%s: impact on the %s correlation, against module properties",
                         GRPLAB[[gr]], PRED[[pv]]),
         subtitle = paste0("y: mean |change in r| over the family's measures at l_max = ", LMF,
                           ". One point per module.\n",
                           "Six points per panel: rho is shown for orientation only and cannot support a claim."),
         x = "property value", y = "mean |change in r|") +
    thm(10) + theme(strip.text.y.left = element_text(angle = 0, hjust = 0))
  ggsave(file.path(outdir, sprintf("properties_%s_%s.png", pv, gr)), p,
         width = 15.0, height = 2.6 * length(unique(long$famlab)) + 3.0,
         dpi = 200, bg = "white")
}

cat("wrote", length(list.files(outdir, "\\.png$")), "figures to", outdir, "\n")
for (f in sort(list.files(outdir, "\\.png$"))) cat("  ", f, "\n")
