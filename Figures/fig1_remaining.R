## fig1_remaining.R --------------------------------------------------------
## The Figure 1 panels that need no new simulation.
##
##   fig1a   influence matrices of the strongest and weakest module (RS, PS)
##   fig1e   frustration vs team strength          (frust0 is already in
##           every finFlagFreq.csv -- no new computation needed)
##   fig1f   hybrid-state frequency F_hy by module
##   fig1j   F_hy vs module team strength
##
## NOT here, because the single-module Ising runs for THIS paper's module
## definitions do not exist on this mac or on Explorer:
##   1c  steady-state heatmap of one module
##   1d  P_s violin over the six modules
## and 1g (P_mean) needs a perturbation simulation that has never been run.
## See FIGURE_PLAN.md "Blocked on the cluster".
##
## Conventions: l_max = 10, Ising for everything, RACIPE for PCA only.
## Run from the repo root:  Rscript Figures/fig1_remaining.R
## -------------------------------------------------------------------------
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(ggrepel)
  library(patchwork); library(funcsKishore)
})
repo <- getwd()
options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))

outdir <- file.path("Figures", "panels")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
LMAX <- 10

MODCOL <- c(AS="#2a78d6", C="#1baf7a", E="#eda100", M="#008300",
            PS="#4a3aa7", RS="#e34948")
NMODCOL <- c(`2`="#86b6ef", `3`="#5598e7", `4`="#2a78d6", `5`="#184f95")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; RED <- "#e34948"
GREY <- "#f0efec"; INK <- "#0b0b0b"; MUTED <- "#52514e"

SRC_ISING <- "Phenotypic scores and frustration: Ising steady states (finFlagFreq, basin-weighted)."
SRC_TOPO  <- "Structural metrics from the influence matrix; no simulation."
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
sv <- function(p, f, w, h) {
  ggsave(file.path(outdir, f), p, width = w, height = h, dpi = 400, bg = "white")
  cat(sprintf("  %-34s %4.1f x %4.1f in\n", f, w, h))
}

N   <- read.csv("Metrics/metrics_network_level.csv", stringsAsFactors = FALSE)
M   <- read.csv("Metrics/metrics_module_level.csv",  stringsAsFactors = FALSE)
PR  <- read.csv("Metrics/module_properties.csv",     stringsAsFactors = FALSE)
n10 <- N[N$l_max == LMAX, ]; m10 <- M[M$l_max == LMAX, ]
pr10 <- PR[PR$l_max == LMAX, ]
MODORD <- pr10$module[order(pr10$T_s_module)]

freq_file <- function(net) {
  for (d in c("Two_module_nets", NETWORK_DIR)) {
    p <- file.path(repo, d, paste0(net, "_finFlagFreq.csv"))
    if (file.exists(p)) return(p)
  }
  NA_character_
}

## =========================================================================
## Fig 1a -- influence matrices of the strongest and weakest module
## RS (T_s 0.92) against PS (0.086). outline.md proposed AS vs PS; RS is the
## sharper contrast and AS is second-strongest, so RS is used and both
## numbers go in the caption.
## =========================================================================
infl_df <- function(m) {
  topo <- file.path(repo, "Single_EMT_MechSens", paste0("Single_", m, ".topo"))
  I  <- influence_full(topo, LMAX)
  tm <- module_teams(repo, m)
  ord <- c(intersect(tm[[1]], rownames(I)), intersect(tm[[2]], rownames(I)))
  I  <- I[ord, ord, drop = FALSE]
  d  <- as.data.frame(as.table(I)); names(d) <- c("from", "to", "v")
  d$from <- factor(d$from, levels = ord); d$to <- factor(d$to, levels = ord)
  d$module <- sprintf("%s   (T_s = %.2f, %d nodes)", m,
                      pr10$T_s_module[pr10$module == m], length(ord))
  d$split <- length(intersect(tm[[1]], ord)) + 0.5
  d
}
## one plot per module: the two modules differ in size (7 vs 15 nodes), so
## they need free scales, which rules out facetting with coord_equal
mk1a <- function(m, show_legend) {
  d <- infl_df(m); sp <- d$split[1]
  ggplot(d, aes(from, to, fill = v)) +
    geom_tile() +
    geom_hline(yintercept = sp, colour = INK, linewidth = 0.6) +
    geom_vline(xintercept = sp, colour = INK, linewidth = 0.6) +
    scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                         limits = c(-1, 1), name = "influence") +
    scale_y_discrete(limits = rev) +
    coord_equal() +
    labs(title = sprintf("%s   (T_s = %.2f)", m, pr10$T_s_module[pr10$module == m]),
         x = "source", y = "target") +
    thm(10) +
    theme(axis.text = element_text(size = rel(0.55)),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          panel.grid.major = element_blank(),
          legend.position = if (show_legend) "right" else "none",
          plot.title = element_text(size = rel(0.95), face = "bold"))
}
p1a <- (mk1a("RS", FALSE) | mk1a("PS", TRUE)) +
  plot_annotation(
    title = "A strong-team and a weak-team module",
    subtitle = "influence matrix, nodes ordered by team; black lines mark the team boundary",
    caption = cap("l_max = 10.", SRC_TOPO,
                  "RS is the strongest module and PS the weakest; AS (0.77) is second-strongest."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(colour = MUTED, size = 9),
                  plot.caption = element_text(colour = MUTED, size = 8, hjust = 0)))
sv(p1a, "fig1a_influence_RS_PS.png", 9.5, 5.0)

## =========================================================================
## Fig 1e -- frustration
## frust0 is already a column of every finFlagFreq.csv, written by the Boolean
## run. Nothing new is computed here: this is the basin-weighted mean over the
## states of each network.
## =========================================================================
frust <- bind_rows(lapply(n10$network, function(net) {
  fp <- freq_file(net); if (is.na(fp)) return(NULL)
  fr <- read.csv(fp, stringsAsFactors = FALSE)
  data.frame(network = net,
             frust_mean = weighted.mean(fr$frust0, fr$Avg0),
             frust_dom  = fr$frust0[which.max(fr$Avg0)],
             n_states   = nrow(fr))
})) %>% left_join(n10[, c("network", "T_network", "n_modules")], by = "network") %>%
  mutate(nmod = factor(n_modules))

p1e <- ggplot(frust, aes(T_network, frust_mean)) +
  geom_smooth(method = "lm", se = TRUE, colour = INK, fill = "grey80",
              linewidth = 0.8, alpha = 0.3) +
  geom_point(aes(fill = nmod), shape = 21, size = 2.6, colour = "white", stroke = 0.7) +
  annotate("text", x = max(frust$T_network), y = max(frust$frust_mean),
           label = sprintf("r = %+.3f", cor(frust$T_network, frust$frust_mean)),
           hjust = 1, vjust = 1, size = 3.2) +
  scale_fill_manual(values = NMODCOL, name = "modules") +
  labs(title = "Stronger teams give less frustrated steady states",
       subtitle = "basin-weighted mean frustration of the Ising attractors, 56 networks",
       x = expression("network team strength  "*T[s]), y = "mean frustration",
       caption = cap("l_max = 10.", SRC_ISING,
                     "Frustration is the frust0 column, written by the Boolean run.")) +
  thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(p1e, "fig1e_frustration.png", 5.6, 4.2)

## =========================================================================
## Fig 1f / 1j -- hybrid states
## A state is hybrid for a module when that module's teams are both partly on,
## i.e. |P_s^state| is far from its committed value of 2. Threshold |P_s| < 1
## (less than halfway committed); the sensitivity over 0.5/1.0/1.5 is written
## to Figures/panels/fig1f_hybrid_sensitivity.csv.
## =========================================================================
THRESH <- c(0.5, 1.0, 1.5); MAIN_T <- 1.0
hyb <- bind_rows(lapply(seq_len(nrow(m10)), function(i) {
  net <- m10$network[i]; mod <- m10$module[i]
  fp <- freq_file(net); if (is.na(fp)) return(NULL)
  ## the state string is ordered by the nodes file sitting next to the
  ## finFlagFreq, so resolve it from there rather than from the topo dir
  nodes <- read_nodes(file.path(dirname(fp), paste0(net, "_nodes.txt")))
  tb <- state_ps_table(fp, nodes, module_teams(repo, mod))
  data.frame(network = net, module = mod,
             thresh = THRESH,
             F_hy = vapply(THRESH, function(th) sum(tb$freq[abs(tb$ps) < th]), numeric(1)))
}))
write.csv(hyb, file.path(outdir, "fig1f_hybrid_sensitivity.csv"), row.names = FALSE)

h <- hyb[hyb$thresh == MAIN_T, ] %>%
  left_join(m10[, c("network", "module", "T_module_intra")],
            by = c("network", "module")) %>%
  mutate(module = factor(module, levels = MODORD))

hs <- hyb %>%
  left_join(m10[, c("network", "module", "T_module_intra")],
            by = c("network", "module")) %>%
  left_join(pr10[, c("module", "n_nodes", "T_s_module")], by = "module") %>%
  mutate(module = factor(module, levels = MODORD),
         thlab  = factor(sprintf("|P_s| < %.1f", thresh)))

## --- 1f: F_hy by module, at all three thresholds -------------------------
## The ranking is not stable across thresholds, which is the finding.
p1f <- ggplot(hs, aes(module, F_hy, fill = module)) +
  geom_violin(colour = NA, alpha = 0.55, width = 0.9, scale = "width") +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white",
               colour = INK, linewidth = 0.35) +
  facet_wrap(~thlab, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = MODCOL, guide = "none") +
  labs(title = "Hybrid-state frequency does not order by team strength",
       subtitle = "modules on x ordered weakest to strongest team; the ranking reshuffles with the threshold",
       x = "module (ordered by team strength)", y = expression(F[hy]),
       caption = cap("l_max = 10.", SRC_ISING,
                     "RS, the *strongest* module, is the most hybrid at every threshold.",
                     "Thresholds are in +/-1 coding (|P_s| = 2 is fully committed); 0.5 / 1.0 / 1.5 here are",
                     "0.25 / 0.5 / 0.75 in the formula sheet's activity-fraction units.",
                     "\nPer-instance values in fig1f_hybrid_sensitivity.csv.")) +
  thm()
sv(p1f, "fig1f_hybrid_by_module.png", 9.0, 4.2)

## --- 1j: what F_hy actually tracks --------------------------------------
agg <- hs %>% group_by(thlab, module) %>%
  summarise(F = mean(F_hy), n_nodes = first(n_nodes),
            Ts = first(T_s_module), .groups = "drop")
lab <- function(v, x, y) sprintf("r = %+.2f", cor(x, y))
ann <- bind_rows(
  agg %>% group_by(thlab) %>% summarise(what = "team strength",
            r = cor(F, Ts), .groups = "drop"),
  agg %>% group_by(thlab) %>% summarise(what = "module size (nodes)",
            r = cor(F, n_nodes), .groups = "drop"))

long <- bind_rows(
  transform(agg, x = Ts,      what = "team strength"),
  transform(agg, x = n_nodes, what = "module size (nodes)"))
long$what <- factor(long$what, levels = c("team strength", "module size (nodes)"))
ann$what  <- factor(ann$what,  levels = levels(long$what))

p1j <- ggplot(long, aes(x, F)) +
  geom_smooth(method = "lm", se = FALSE, colour = INK, linewidth = 0.7) +
  geom_point(aes(fill = module), shape = 21, size = 3, colour = "white", stroke = 0.7) +
  geom_text_repel(aes(label = module, colour = module), size = 2.6,
                  fontface = "bold", seed = 1, box.padding = 0.3,
                  min.segment.length = 0.3, segment.size = 0.2, show.legend = FALSE) +
  geom_text(data = ann, aes(x = Inf, y = Inf, label = sprintf("r = %+.2f", r)),
            hjust = 1.12, vjust = 1.6, size = 3.0, colour = MUTED, inherit.aes = FALSE) +
  facet_grid(thlab ~ what, scales = "free") +
  scale_fill_manual(values = MODCOL, guide = "none") +
  scale_colour_manual(values = MODCOL, guide = "none") +
  labs(title = expression(F[hy]*" tracks module size, not team strength"),
       subtitle = "six module means, at each hybrid threshold; smaller modules look more hybrid",
       x = NULL, y = expression("mean "*F[hy]),
       caption = cap("l_max = 10.", SRC_ISING,
                     "A fixed |P_s| cut is not size-invariant: a 6-node module quantises P_s far more coarsely than a 17-node one,",
                     "\nso it lands in the hybrid band more often. The definition needs revising before F_hy carries any weight.")) +
  thm() + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))
sv(p1j, "fig1j_hybrid_vs_Ts.png", 8.0, 7.0)

cat(sprintf("\n  r(T_network, frustration)      = %+.3f\n", cor(frust$T_network, frust$frust_mean)))
cat("  F_hy, across the six module means:\n")
print(as.data.frame(agg %>% group_by(thlab) %>%
  summarise(r_team_strength = round(cor(F, Ts), 3),
            r_module_size   = round(cor(F, n_nodes), 3), .groups = "drop")))
cat("\ndone.\n")
