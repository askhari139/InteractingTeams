## fig1_build.R -----------------------------------------------------------
## Builds every panel of the REVISED Figure 1 and the compiled figure, to
## the layout in Figures/correlation_heatmaps/figure_outline_new.md:
##
##   Row 1   A                RS: influence matrix + RACIPE states
##   Row 2   B                PS: influence matrix + RACIPE states
##   Row 3   C | D
##   Row 4   E  |  F over G
##
## Figure 1 leads with RACIPE (directive 2026-09-24). No phenotypic score
## appears, and panel D carries NO module-derived predictor.
##
## Conventions set 2026-09-24:
##   - every legend sits on TOP of its panel
##   - the COMPILED figure carries no title, subtitle or caption text; panel
##     identity comes from cowplot's A-G labels. The individually saved
##     panels keep their titles and captions.
##   - compilation uses cowplot::plot_grid, not patchwork
##
## Output: Figures/figure1_new/   Run from the repo root.
## -------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(dplyr); library(tidyr)
                  library(cowplot); library(ggrepel); library(funcsKishore)})
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))
outdir <- file.path("Figures", "figure1_new"); dir.create(outdir, FALSE, TRUE)
zdir <- "RACIPE_zscores"; LMAX <- 10

MODCOL <- c(AS = "#2a78d6", C = "#1baf7a", E = "#eda100", M = "#008300",
            PS = "#4a3aa7", RS = "#e34948")
## Module count in F is drawn as a DISCRETE category, not an ordinal ramp:
## the first five slots of the reference categorical order, which passes every
## adjacent-pair gate (violins are dodged, so adjacent pairs are what matter).
NMODCOL <- c("1" = "#2a78d6", "2" = "#eb6834", "3" = "#1baf7a",
             "4" = "#eda100", "5" = "#e87ba4")
## D's four series: validated all-pairs in light mode (worst CVD dE 9.2,
## worst normal-vision dE 16.3). The aqua sits under 3:1 on the surface, so
## every series is direct-labelled as well as coloured.
DCOL <- c("#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7")
BLUE <- "#2a78d6"; RED <- "#e34948"; GREY <- "#f0efec"
INK <- "#0b0b0b"; MUTED <- "#52514e"

cap <- function(..., width = 115)
  paste(strwrap(paste(Filter(nzchar, c(...)), collapse = " "), width = width),
        collapse = "\n")

## legends on top, laid out horizontally
thm <- function(base = 11) theme_Publication(base_size = base) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        plot.title    = element_text(hjust = 0, size = rel(1.0), face = "bold"),
        plot.subtitle = element_text(hjust = 0, colour = MUTED, size = rel(0.72)),
        plot.caption  = element_text(hjust = 0, colour = MUTED, size = rel(0.62)),
        axis.title    = element_text(size = rel(0.85)),
        legend.position  = "top",
        legend.direction = "horizontal",
        legend.box       = "horizontal",
        legend.justification = "left",
        legend.title  = element_text(size = rel(0.72)),
        legend.text   = element_text(size = rel(0.7)),
        legend.key.size = unit(0.32, "cm"),
        strip.text    = element_text(size = rel(0.8)))

thm_tile <- function(base = 11) thm(base) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank())

## a horizontal colourbar, sized for a legend sitting above the panel
hbar <- function(title) guide_colourbar(
  title = title, title.position = "left", title.vjust = 0.9,
  barwidth = unit(3.0, "cm"), barheight = unit(0.28, "cm"),
  ticks.colour = "white")

## strip every piece of text for the compiled figure
bare <- function(p) p + labs(title = NULL, subtitle = NULL, caption = NULL)

sv <- function(p, f, w, h) {
  ggsave(file.path(outdir, f), p, width = w, height = h, dpi = 400,
         bg = "white", limitsize = FALSE)
  cat(sprintf("  %-32s %5.1f x %5.1f in\n", f, w, h))
}

rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
N <- bind_rows(rd("Metrics/metrics_network_level.csv"),
               rd("Metrics/metrics_network_level_singles.csv"))
M <- bind_rows(rd("Metrics/metrics_module_level.csv"),
               rd("Metrics/metrics_module_level_singles.csv"))
N <- N[N$l_max == LMAX, ]; M <- M[M$l_max == LMAX, ]
PR <- rd("Metrics/module_properties.csv"); pr10 <- PR[PR$l_max == LMAX, ]
MODORD <- pr10$module[order(pr10$T_s_module)]

team_order <- function(mods, pool) {
  ord <- c(); bteam <- c(); bmod <- c()
  for (m in mods) {
    tm <- module_teams(repo, m)
    t1 <- intersect(tm[[1]], pool); t2 <- intersect(tm[[2]], pool)
    st <- length(ord); ord <- c(ord, t1, t2)
    bteam <- c(bteam, st + length(t1) + 0.5)
    bmod  <- c(bmod, length(ord) + 0.5)
  }
  list(order = ord, team = bteam, module = head(bmod, -1))
}

infl_panel <- function(net, mods, show_leg = FALSE, ttl = NULL) {
  tp <- find_net_file(repo, net, "topo")
  old <- setwd(dirname(tp)); I <- influence_full(basename(tp), LMAX); setwd(old)
  to <- team_order(mods, rownames(I)); n <- length(to$order)
  d <- as.data.frame(as.table(I[to$order, to$order, drop = FALSE]))
  names(d) <- c("from", "to", "v")
  d$from <- factor(d$from, levels = to$order)
  d$to   <- factor(d$to,   levels = to$order)
  p <- ggplot(d, aes(from, to, fill = v)) + geom_tile() +
    geom_vline(xintercept = to$team, colour = INK, linewidth = 0.3, linetype = "22") +
    geom_hline(yintercept = n - to$team + 1, colour = INK, linewidth = 0.3, linetype = "22")
  if (length(to$module)) {
    p <- p + geom_vline(xintercept = to$module, colour = INK, linewidth = 0.8) +
      geom_hline(yintercept = n - to$module + 1, colour = INK, linewidth = 0.8)
  }
  ## no coord_equal: with 7 nodes (RS) a square matrix leaves most of the
  ## row empty in the compiled figure. Cells stretch to fill instead.
  p + scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                           limits = c(-1, 1), guide = hbar("influence")) +
    scale_y_discrete(limits = rev) +
    labs(title = ttl, x = "source", y = "target") +
    thm_tile(10) +
    theme(axis.text = element_text(size = rel(0.55)),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = if (show_leg) "top" else "none",
          plot.title = element_text(size = rel(0.95), face = "bold"))
}

racipe_panel <- function(net, mods, show_leg = FALSE, row_cap = 1200, ttl = NULL) {
  d <- read.csv(file.path(zdir, paste0(net, "_z.csv")), check.names = FALSE,
                comment.char = "#")
  genes <- setdiff(names(d), "basin"); to <- team_order(mods, genes)
  Z <- as.matrix(d[, to$order, drop = FALSE])
  t1 <- intersect(module_teams(repo, mods[1])[[1]], to$order)
  t2 <- intersect(module_teams(repo, mods[1])[[2]], to$order)
  keep <- order(rowMeans(Z[, t1, drop = FALSE]) - rowMeans(Z[, t2, drop = FALSE]))
  if (length(keep) > row_cap)
    keep <- keep[round(seq(1, length(keep), length.out = row_cap))]
  Z <- Z[keep, , drop = FALSE]
  rownames(Z) <- sprintf("r%05d", seq_len(nrow(Z)))
  long <- as.data.frame(as.table(Z), stringsAsFactors = FALSE)
  names(long) <- c("row", "gene", "z")
  long$row  <- factor(long$row,  levels = rownames(Z))
  long$gene <- factor(long$gene, levels = to$order)
  p <- ggplot(long, aes(gene, row, fill = z)) + geom_raster() +
    geom_vline(xintercept = to$team, colour = INK, linewidth = 0.3, linetype = "22")
  if (length(to$module))
    p <- p + geom_vline(xintercept = to$module, colour = INK, linewidth = 0.8)
  p + scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                           limits = c(-2.5, 2.5), oob = scales::squish,
                           guide = hbar("z-score")) +
    labs(title = ttl, x = "node (ordered by team)",
         y = sprintf("RACIPE solutions (%d shown)", nrow(Z))) +
    thm_tile(10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(0.55)),
          axis.text.y = element_blank(),
          axis.title = element_text(size = rel(0.7)),
          legend.position = if (show_leg) "top" else "none",
          plot.title = element_text(size = rel(0.95), face = "bold"))
}

## ===== A, B : one module per row, influence matrix + RACIPE states =======
A_infl <- infl_panel("Comb_RS_RS", "RS", TRUE,
  sprintf("A   RS influence matrix  (T_s = %.2f)", pr10$T_s_module[pr10$module == "RS"]))
A_st   <- racipe_panel("Comb_RS_RS", "RS", TRUE, 1200, "A   RS RACIPE state space")
B_infl <- infl_panel("Comb_PS_PS", "PS", TRUE,
  sprintf("B   PS influence matrix  (T_s = %.2f)", pr10$T_s_module[pr10$module == "PS"]))
B_st   <- racipe_panel("Comb_PS_PS", "PS", TRUE, 1200, "B   PS RACIPE state space")
sv(plot_grid(A_infl, A_st, nrow = 1, rel_widths = c(1, 1.25)), "fig1A_RS.png", 11.0, 4.8)
sv(plot_grid(B_infl, B_st, nrow = 1, rel_widths = c(1, 1.25)), "fig1B_PS.png", 11.0, 4.8)

## ===== C : heatmap, x = module, y = metric, fill = value =================
CP <- c(T_s_module = "team strength", density_within = "density",
        impurity_across = "impurity")
dC <- pr10 %>% select(module, all_of(names(CP))) %>%
  pivot_longer(-module, names_to = "prop", values_to = "raw") %>%
  group_by(prop) %>% mutate(z = as.numeric(scale(raw))) %>% ungroup() %>%
  mutate(module = factor(module, levels = MODORD),
         prop = factor(CP[prop], levels = CP))
C <- ggplot(dC, aes(module, prop, fill = z)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%.2f", raw), colour = abs(z) > 1.1),
            size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = INK)) +
  scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                       guide = hbar("z within metric")) +
  scale_y_discrete(limits = rev(levels(dC$prop))) +
  labs(title = "C   Module properties",
       x = "module (ordered by team strength)", y = NULL,
       caption = cap("l_max = 10; influence-matrix metrics, no simulation.",
                     "Cells give the raw value; fill is scaled within each metric row.",
                     "Impurity is the cross-module value: within-module impurity is 0",
                     "for every module but PS (0.143).")) +
  thm_tile()
sv(C, "fig1C_module_properties.png", 7.4, 3.2)

## ===== D : PC1 variance and PCA entropy against team strength ============
## Both outcomes are already on a 0-1 scale, so they share the y axis at
## their RAW values -- no z-scoring, and no second axis.
## Colour encodes the y-axis metric, not the predictor.
DOUT <- c(PC1_var = "PC1 variance", PCA_entropy_norm = "PCA entropy (norm)")
dD <- N %>% select(network, n_modules, T_network, all_of(names(DOUT))) %>%
  pivot_longer(all_of(names(DOUT)), names_to = "metric", values_to = "v") %>%
  mutate(metric = factor(DOUT[metric], levels = DOUT))
rl <- dD %>% group_by(metric) %>%
  summarise(r = cor(T_network, v), x = max(T_network), y = max(v), .groups = "drop") %>%
  mutate(lab = sprintf("%s  r = %+.2f", metric, r))
D <- ggplot(dD, aes(T_network, v, colour = metric)) +
  geom_point(size = 1.6, alpha = 0.6, stroke = 0) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  geom_text_repel(data = rl, aes(x, y, label = lab), size = 2.9, fontface = "bold",
                  seed = 1, direction = "y", hjust = 1, box.padding = 0.5,
                  segment.size = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = DCOL[1:2], name = NULL) +
  labs(title = "D   Stronger teams concentrate the RACIPE state space",
       subtitle = sprintf("one point per network (n = %d); colour is the y-axis metric", nrow(N)),
       x = expression("network team strength  " * T[s]),
       y = "metric value",
       caption = cap("PCA metrics from RACIPE continuous sampling; team strength is an",
                     "influence-matrix metric at l_max = 10. Both metrics are already on a",
                     "0-1 scale, so they share the axis at their raw values with no rescaling.",
                     "The two run in opposite directions: concentration up, diversity down.")) +
  thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(D, "fig1D_pc_vs_Ts.png", 7.4, 4.6)

## ===== E : blank third on the right, for manual annotation ===============
E <- infl_panel("Comb_PS_RS", c("PS", "RS"), TRUE,
                "E   Two coupled modules: PS + RS") +
  labs(caption = cap("l_max = 10. Solid lines: module boundary.",
                     "Dashed: team split within a module.",
                     "Off-diagonal blocks are the inter-module coupling."))
sv(plot_grid(E, NULL, nrow = 1, rel_widths = c(2, 1)),
   "fig1E_influence_PS_RS.png", 10.0, 5.6)

## ===== F : intra-module T_s by module count ==============================
dF <- M %>% mutate(module = factor(module, levels = MODORD),
                   nmod = factor(n_modules))
F1 <- ggplot(dF, aes(module, T_module_intra, fill = nmod)) +
  geom_violin(position = position_dodge(width = 0.8), colour = NA,
              alpha = 0.75, scale = "width", width = 0.75) +
  ## The median marker carries the fill too. n = 1 is a single value per
  ## module, so geom_violin drops that group entirely -- without a filled
  ## marker its legend key would render blank.
  stat_summary(fun = median, geom = "point",
               position = position_dodge(width = 0.8),
               shape = 21, size = 1.9, colour = INK, stroke = 0.4) +
  scale_fill_manual(values = NMODCOL, name = "modules", drop = FALSE) +
  labs(title = "F   A module's team strength shifts with how many modules it sits in",
       x = "module (ordered by team strength in isolation)",
       y = expression("intra-module  " * T[s]^{mod}),
       caption = cap("l_max = 10. n = 1 is the module alone, so it has a single value.",
                     "The module's own teams are evaluated on the whole network's",
                     "influence matrix, which is why the value moves at all.")) +
  thm()
sv(F1, "fig1F_intra_Ts_by_nmodules.png", 8.4, 4.2)

## ===== G : PC metrics against module count ===============================
PCM <- c(PC1_var = "PC1 variance", NumPC_90 = "PCs to 90%",
         PCA_entropy_norm = "PCA entropy (norm)")
dG <- N %>% select(network, n_modules, all_of(names(PCM))) %>%
  pivot_longer(-c(network, n_modules), names_to = "pc", values_to = "v") %>%
  group_by(pc) %>% mutate(z = as.numeric(scale(v))) %>% ungroup() %>%
  mutate(pc = factor(PCM[pc], levels = PCM), nmod = factor(n_modules))
gm <- dG %>% group_by(pc, nmod) %>% summarise(z = mean(z), .groups = "drop")
G <- ggplot(dG, aes(nmod, z, colour = pc, fill = pc)) +
  geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.35) +
  geom_violin(position = position_dodge(width = 0.8), colour = NA, alpha = 0.3,
              scale = "width", width = 0.75) +
  geom_line(data = gm, aes(group = pc), position = position_dodge(width = 0.8),
            linewidth = 0.8) +
  geom_point(data = gm, position = position_dodge(width = 0.8), size = 2) +
  scale_colour_manual(values = DCOL[1:3], name = NULL) +
  scale_fill_manual(values = DCOL[1:3], name = NULL) +
  labs(title = "G   The RACIPE state space diversifies as modules are added",
       x = "modules in the network", y = "z-score (within metric)",
       caption = cap("Each PC metric z-scored so the three share one axis;",
                     "lines join the per-group means.")) +
  thm()
sv(G, "fig1G_pc_vs_nmodules.png", 8.4, 4.2)

## ===== supplementary: every module in isolation ==========================
sp <- lapply(seq_along(MODORD), function(i) {
  m <- MODORD[i]
  racipe_panel(sprintf("Comb_%s_%s", m, m), m, show_leg = TRUE, row_cap = 900,
               sprintf("%s  (T_s = %.2f)", m, pr10$T_s_module[pr10$module == m])) +
    labs(x = NULL, y = if (i %% 3 == 1) "RACIPE solutions" else NULL)
})
sv(plot_grid(plotlist = sp, nrow = 2), "figS_racipe_all_modules.png", 15.0, 9.0)

## ===== compiled figure ===================================================
## No title, subtitle or caption anywhere; panel identity is carried by the
## cowplot A-G labels. Assembled with cowplot::plot_grid.
LS <- 17
row1 <- plot_grid(bare(A_infl), bare(A_st), nrow = 1, rel_widths = c(1, 1.25),
                  labels = c("A", ""), label_size = LS)
row2 <- plot_grid(bare(B_infl), bare(B_st), nrow = 1, rel_widths = c(1, 1.25),
                  labels = c("B", ""), label_size = LS)
row3 <- plot_grid(bare(C), bare(D), nrow = 1, rel_widths = c(1, 1.1),
                  labels = c("C", "D"), label_size = LS)
fg   <- plot_grid(bare(F1), bare(G), ncol = 1, labels = c("F", "G"), label_size = LS)
row4 <- plot_grid(plot_grid(bare(E), NULL, nrow = 1, rel_widths = c(2, 1)), fg,
                  nrow = 1, rel_widths = c(1, 1), labels = c("E", ""), label_size = LS)
FIG <- plot_grid(row1, row2, row3, row4, ncol = 1,
                 rel_heights = c(1.0, 1.0, 0.95, 1.7))
sv(FIG, "FIG1_compiled.png", 17.5, 23.0)
cat("\ndone.\n")
