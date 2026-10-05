## fig1_singles.R ---------------------------------------------------------
## Figure 1 panels that need the single-module Ising runs.
##
## The single-module simulations are the SELF-PAIRS in Two_module_nets:
## Comb_AS_AS, Comb_C_C, ... Each is byte-identical in topology and node order
## to Single_<M>, verified 2026-09-22 -- so no new simulation is needed.
## compute_metrics.R deliberately skips them ("drop self-pairs"), which is why
## they were missed the first time round.
##
##   fig1c   steady-state heatmap (states x nodes, team-ordered) + P_s strip
##   fig1d   per-state P_s by module, basin-weighted, ordered by team strength
##
## Conventions: Ising for everything; l_max = 10 only for the T_s ordering.
## Run from the repo root:  Rscript Figures/fig1_singles.R
## -------------------------------------------------------------------------
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(ggrepel); library(funcsKishore)
})
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))
outdir <- file.path("Figures", "panels"); dir.create(outdir, FALSE, TRUE)

MODCOL <- c(AS="#2a78d6", C="#1baf7a", E="#eda100", M="#008300",
            PS="#4a3aa7", RS="#e34948")
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; RED <- "#e34948"
GREY <- "#f0efec"; INK <- "#0b0b0b"; MUTED <- "#52514e"
SRC <- "Single-module Ising runs: the Comb_<M>_<M> self-pairs in Two_module_nets (topology identical to Single_<M>)."
cap <- function(..., width = 110)
  paste(strwrap(paste(c(...), collapse = " "), width = width), collapse = "\n")

thm <- function(base = 11) theme_Publication(base_size = base) +
  theme(panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        plot.title    = element_text(hjust = 0, size = rel(1.0), face = "bold"),
        plot.subtitle = element_text(hjust = 0, colour = MUTED, size = rel(0.72)),
        plot.caption  = element_text(hjust = 0, colour = MUTED, size = rel(0.62)),
        axis.title    = element_text(size = rel(0.92)),
        legend.text   = element_text(size = rel(0.72)),
        legend.title  = element_text(size = rel(0.75)),
        legend.key.size = unit(0.35, "cm"))
sv <- function(p, f, w, h) {
  ggsave(file.path(outdir, f), p, width = w, height = h, dpi = 400, bg = "white")
  cat(sprintf("  %-34s %4.1f x %4.1f in\n", f, w, h))
}

MODS <- c("AS", "C", "E", "M", "PS", "RS")
PR   <- read.csv("Metrics/module_properties.csv"); pr10 <- PR[PR$l_max == 10, ]
MODORD <- pr10$module[order(pr10$T_s_module)]
selfnet <- function(m) sprintf("Comb_%s_%s", m, m)

## per-state table for one module's own simulation
single_states <- function(m) {
  net <- selfnet(m)
  fp  <- file.path("Two_module_nets", paste0(net, "_finFlagFreq.csv"))
  nd  <- readLines(file.path("Two_module_nets", paste0(net, "_nodes.txt")))
  tm  <- module_teams(repo, m)
  fr  <- read.csv(fp, stringsAsFactors = FALSE)
  ps  <- vapply(fr$states, function(s)
           pheno_score_state_based(decode_state(s, nd), tm), numeric(1))
  list(freq = fr$Avg0, ps = as.numeric(ps), states = fr$states,
       nodes = nd, teams = tm, frust = fr$frust0)
}
S <- setNames(lapply(MODS, single_states), MODS)
cat("states per module: ",
    paste(sprintf("%s=%d", MODS, vapply(S, function(z) length(z$freq), 1L)),
          collapse = "  "), "\n\n")

## =========================================================================
## Fig 1c -- steady states of a strong-team and a weak-team module
## Same pairing as 1a (RS vs PS) so the two panels tell one story.
## =========================================================================
one_1c <- function(m, show_leg) {
  z <- S[[m]]
  ord <- order(-z$freq)                       # states by basin size
  tnodes <- c(intersect(z$teams[[1]], z$nodes), intersect(z$teams[[2]], z$nodes))
  mat <- t(vapply(z$states[ord], function(s) decode_state(s, z$nodes)[tnodes],
                  numeric(length(tnodes))))
  d <- as.data.frame(as.table(as.matrix(mat)))
  names(d) <- c("state", "node", "v")
  d$state <- factor(d$state, levels = rev(z$states[ord]))
  d$node  <- factor(d$node,  levels = tnodes)
  nb <- length(intersect(z$teams[[1]], z$nodes)) + 0.5

  hm <- ggplot(d, aes(node, state, fill = factor(v))) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_vline(xintercept = nb, colour = INK, linewidth = 0.7) +
    scale_fill_manual(values = c(`-1` = "#cde2fb", `1` = "#184f95"),
                      labels = c("off", "on"), name = "node state") +
    labs(title = sprintf("%s   (T_s = %.2f, %d states)", m,
                         pr10$T_s_module[pr10$module == m], length(z$freq)),
         x = "node (ordered by team)", y = "steady state (by basin size)") +
    thm(10) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                               size = rel(0.6)),
                    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    panel.grid.major = element_blank(),
                    legend.position = if (show_leg) "right" else "none",
                    plot.title = element_text(size = rel(0.95), face = "bold"))
  st <- data.frame(state = factor(z$states[ord], levels = rev(z$states[ord])),
                   ps = z$ps[ord], freq = z$freq[ord])
  strip <- ggplot(st, aes(1, state, fill = ps)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0,
                         limits = c(-2, 2), name = expression(P[s])) +
    labs(x = NULL, y = NULL, title = " ") +
    thm(10) + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                    panel.grid.major = element_blank(),
                    legend.position = if (show_leg) "right" else "none")
  bas <- ggplot(st, aes(freq, state)) +
    geom_col(fill = "grey55", width = 0.75) +
    scale_x_continuous(breaks = scales::breaks_pretty(2)) +
    labs(x = "basin", y = NULL, title = " ") +
    thm(10) + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    panel.grid.major.y = element_blank())
  hm + strip + bas + plot_layout(widths = c(length(tnodes), 1.1, 3))
}
p1c <- (one_1c("RS", FALSE) / one_1c("PS", TRUE)) +
  plot_annotation(
    title = "Both modules have a dominant mirror pair -- only the strong-team one is fully committed",
    subtitle = "left: steady states, nodes ordered by team (black line = team boundary)  |  middle: that state's phenotypic score  |  right: basin size",
    caption = cap(SRC,
      "P_s = mean(team 1) - mean(team 2) in +/-1 coding, range [-2, 2]. This is exactly 2 x (E(T1) - E(T2)) in the formula sheet's",
      "activity-fraction units, where it would range [-1, 1]; halve these numbers to read them there.",
      "RS: 6 states, top pair holds 78% of the basin at |P_s| = 2.00.",
      "PS: 16 states, top pair holds 59% at |P_s| = 1.83, and only 14% of its basin is fully committed",
      "(against 58-100% for every other module)."),
    theme = theme(plot.title = element_text(face = "bold", size = 13),
                  plot.subtitle = element_text(colour = MUTED, size = 8.5),
                  plot.caption = element_text(colour = MUTED, size = 8, hjust = 0)))
sv(p1c, "fig1c_steady_states.png", 11.0, 7.5)

## =========================================================================
## Fig 1d -- per-state P_s, by module, ordered by team strength
## Each module has only 2-16 states, so a violin would be a lie: plot the
## states themselves, sized by basin, with the basin-weighted mean marked.
## =========================================================================
d1d <- bind_rows(lapply(MODS, function(m) data.frame(
  module = m, ps = S[[m]]$ps, freq = S[[m]]$freq))) %>%
  mutate(module = factor(module, levels = MODORD))
mu <- d1d %>% group_by(module) %>%
  summarise(wm = weighted.mean(abs(ps), freq),
            n = n(), .groups = "drop")

p1d <- ggplot(d1d, aes(module, ps)) +
  geom_hline(yintercept = c(-2, 0, 2), colour = "grey85", linewidth = 0.35) +
  geom_point(aes(size = freq, colour = module), alpha = 0.75) +
  geom_point(data = mu, aes(module, wm), shape = 95, size = 11, colour = INK) +
  geom_text(data = mu, aes(module, 2.42, label = sprintf("%d states", n)),
            size = 2.7, colour = MUTED) +
  scale_colour_manual(values = MODCOL, guide = "none") +
  scale_size_area(max_size = 7, name = "basin size",
                  breaks = c(0.05, 0.25, 0.5)) +
  scale_y_continuous(limits = c(-2.15, 2.6), breaks = seq(-2, 2, 1)) +
  labs(title = "In isolation every module commits, whatever its team strength",
       subtitle = "every steady state of each module's own Ising run; black bar is the basin-weighted mean |P_s|",
       x = "module (ordered by team strength)",
       y = expression("state phenotypic score  "*P[s]),
       caption = cap(SRC,
         "+/-2 means fully committed to one team. The 6 modules have 2-16 states each, so the states are drawn",
         "individually rather than as a violin. +/-1 coding, so |P_s| = 2 is fully committed (= 1 in formula-sheet units).",
         "Mean |P_s| spans only 1.72-2.00 and does not order by team strength",
         "(r = -0.10, n = 6) -- the team-strength dependence appears only once modules are combined (panel k).")) +
  thm()
sv(p1d, "fig1d_state_Ps_by_module.png", 7.0, 4.6)

## =========================================================================
## Fig 1k -- [dev] not in outline.md. Isolated vs embedded commitment.
## 1d shows every module commits on its own. The question the paper asks is
## what combination costs, so put the two side by side.
## =========================================================================
M10 <- read.csv("Metrics/metrics_module_level.csv")
M10 <- M10[M10$l_max == 10, ]
ff <- function(net) {
  for (d in c("Two_module_nets", NETWORK_DIR)) {
    q <- file.path(repo, d, paste0(net, "_finFlagFreq.csv"))
    if (file.exists(q)) return(q)
  }
  NA_character_
}
wmabs <- function(net, m) {
  fp <- ff(net); nd <- readLines(file.path(dirname(fp), paste0(net, "_nodes.txt")))
  tb <- state_ps_table(fp, nd, module_teams(repo, m))
  weighted.mean(abs(tb$ps), tb$freq)
}
emb <- bind_rows(lapply(MODS, function(m) {
  nets <- M10$network[M10$module == m]
  data.frame(module = m, network = nets,
             wm = vapply(nets, function(n) wmabs(n, m), numeric(1)))
}))
iso <- data.frame(module = MODS,
                  iso = vapply(MODS, function(m) wmabs(selfnet(m), m), numeric(1)),
                  T_s = pr10$T_s_module[match(MODS, pr10$module)])
embm <- emb %>% group_by(module) %>% summarise(embm = mean(wm), .groups = "drop")
cmp <- left_join(iso, embm, by = "module") %>% mutate(drop = iso - embm)
cmp$module <- factor(cmp$module, levels = MODORD)
emb$module <- factor(emb$module, levels = MODORD)

pk1 <- ggplot(emb, aes(module, wm)) +
  geom_violin(aes(fill = module), colour = NA, alpha = 0.45, scale = "width", width = 0.85) +
  geom_point(data = cmp, aes(module, iso), shape = 18, size = 3.6, colour = INK) +
  geom_point(data = cmp, aes(module, embm), shape = 21, size = 2.6,
             fill = "white", colour = INK, stroke = 0.7) +
  scale_fill_manual(values = MODCOL, guide = "none") +
  labs(title = "What embedding costs each module",
       subtitle = "violin: the 30 networks containing the module  |  diamond: the module alone  |  open circle: mean when embedded",
       x = "module (ordered by team strength)",
       y = expression("basin-weighted mean  "*group("|",P[s],"|"))) +
  thm(10)

pk2 <- ggplot(cmp, aes(T_s, drop)) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_smooth(method = "lm", se = FALSE, colour = INK, linewidth = 0.7) +
  geom_point(aes(fill = module), shape = 21, size = 3.4, colour = "white", stroke = 0.8) +
  geom_text_repel(aes(label = module, colour = module), size = 2.9, fontface = "bold",
                  seed = 1, box.padding = 0.35, min.segment.length = 0.3,
                  segment.size = 0.25, show.legend = FALSE) +
  annotate("text", x = Inf, y = Inf, hjust = 1.15, vjust = 1.8, size = 3.0, colour = MUTED,
           label = sprintf("r = %+.2f", cor(cmp$T_s, cmp$drop))) +
  scale_fill_manual(values = MODCOL, guide = "none") +
  scale_colour_manual(values = MODCOL, guide = "none") +
  labs(title = "Weak teams lose the most",
       subtitle = "loss of commitment on embedding, against team strength",
       x = expression("module team strength  "*T[s]^{mod}),
       y = expression("isolated "-" embedded")) +
  thm(10) + theme(panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3))

p1k <- (pk1 | pk2) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Team strength does not set commitment in isolation -- it sets how much combination costs",
    caption = cap(SRC,
      "Only 6 modules, so r = -0.51 is suggestive, not significant. M is the exception: it gains commitment when embedded."),
    theme = theme(plot.title = element_text(face = "bold", size = 12.5),
                  plot.caption = element_text(colour = MUTED, size = 8, hjust = 0)))
sv(p1k, "fig1k_isolated_vs_embedded.png", 11.0, 4.6)

cat("
  isolated vs embedded commitment:
")
print(as.data.frame(cmp[order(cmp$T_s), ]), row.names = FALSE, digits = 3)
cat(sprintf("  r(T_s, isolated) = %+.3f   r(T_s, embedded) = %+.3f   r(T_s, drop) = %+.3f\n",
    cor(cmp$T_s, cmp$iso), cor(cmp$T_s, cmp$embm), cor(cmp$T_s, cmp$drop)))

cat("\n  basin-weighted mean |P_s| by module (weakest to strongest team):\n")
print(as.data.frame(mu))
cat(sprintf("\n  r(module T_s, weighted mean |P_s|) = %+.3f\n",
    cor(pr10$T_s_module[match(mu$module, pr10$module)], mu$wm)))
cat("\ndone.\n")
