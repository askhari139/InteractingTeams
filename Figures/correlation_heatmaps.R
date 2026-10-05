## correlation_heatmaps.R -------------------------------------------------
## Input x output correlation matrices, faceted by module count (1-5), at
## l_max = 10. Classified by OUTPUT FAMILY so each matrix has one grain and
## one readout:
##
##   1  corr_pc_lmax10.png            PC + variant-independent metrics
##   2  corr_module_infl_lmax10.png   module P_s, influence-weighted
##   3  corr_module_state_lmax10.png  module P_s, state-based
##   4  corr_network_infl_lmax10.png  network P_s, influence-weighted
##   5  corr_network_state_lmax10.png network P_s, state-based
##   6  corr_modulegrain_lmax10.png   TRUE module grain -- see below
##
## Matrices 1-5 are at NETWORK grain: one row per network. A network has 1-5
## modules but only one value of each structural input, so module P_s is
## rolled up. Both roll-ups are reported:
##   MEAN over the network's modules -- reproduces the published *_module
##        rows exactly (to 1.9e-5, the rounding in those files)
##   VARIANCE over the network's modules -- asks whether structure predicts
##        how UNEVENLY the modules score, which the mean cannot see.
##        Undefined (NA) at n = 1, where there is one module.
##
## Matrix 6 avoids the roll-up entirely: one row per (network, module), with
## module-grain inputs (that module's own intra-T_s and its in/out/total
## inter-module coupling) against that module's own P_s. Sample sizes are
## larger than the network-grain matrices: 6/30/60/60/30.
##
## Inter-module inputs use the TEAM-AWARE metric (T_inter_blocks, see
## metrics_lib.R): the old whole-block T_inter cancels team-aligned coupling
## to ~0. The old metric is kept alongside as a labelled row for continuity.
##
## Run from the repo root:  Rscript Figures/correlation_heatmaps.R
## -------------------------------------------------------------------------
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(funcsKishore)})
repo <- getwd(); options(interactingteams.repo = repo)
source(file.path(repo, "Metrics", "metrics_lib.R"))
outdir <- file.path("Figures", "correlation_heatmaps")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
unlink(Sys.glob(file.path(outdir, "corr_heatmap_*")))   # retire the old 3-file layout
LMAX <- 10

BLUE <- "#2a78d6"; RED <- "#e34948"; GREY <- "#f0efec"
INK <- "#0b0b0b"; MUTED <- "#52514e"
cap <- function(..., width = 150)
  paste(strwrap(paste(Filter(nzchar, c(...)), collapse = " "), width = width), collapse = "\n")
thm <- function(base = 11) theme_Publication(base_size = base) +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.ticks = element_blank(),
        plot.title    = element_text(hjust = 0, size = rel(1.0), face = "bold"),
        plot.subtitle = element_text(hjust = 0, colour = MUTED, size = rel(0.72)),
        plot.caption  = element_text(hjust = 0, colour = MUTED, size = rel(0.62)),
        strip.text    = element_text(size = rel(0.8)),
        legend.key.size = unit(0.35, "cm"))

## ---- data ---------------------------------------------------------------
rd <- function(f) if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
N  <- bind_rows(rd("Metrics/metrics_network_level.csv"), rd("Metrics/metrics_network_level_singles.csv"))
M  <- bind_rows(rd("Metrics/metrics_module_level.csv"),  rd("Metrics/metrics_module_level_singles.csv"))
TA <- rd("Metrics/metrics_inter_teamaware.csv")
TN <- rd("Metrics/metrics_inter_teamaware_network.csv")
N <- N[N$l_max == LMAX, ]; M <- M[M$l_max == LMAX, ]
TA <- TA[TA$l_max == LMAX, ]; TN <- TN[TN$l_max == LMAX, ]
stopifnot(nrow(N) == 62)

## frustration + state count, from the frust0 column already in finFlagFreq
ff <- function(net) { for (d in c("Two_module_nets", NETWORK_DIR)) {
  p <- file.path(repo, d, paste0(net, "_finFlagFreq.csv")); if (file.exists(p)) return(p) }; NA_character_ }
FR <- bind_rows(lapply(N$network, function(net) {
  fr <- read.csv(ff(net), stringsAsFactors = FALSE)
  data.frame(network = net, frustration = weighted.mean(fr$frust0, fr$Avg0), n_states = nrow(fr)) }))

## module P_s rolled up two ways
roll <- M %>% group_by(network) %>%
  summarise(across(starts_with("P_s_"), list(mean = ~mean(.x), var = ~stats::var(.x)),
                   .names = "{.col}__{.fn}"), .groups = "drop")

D <- N %>% left_join(roll, by = "network") %>% left_join(FR, by = "network") %>%
  left_join(TN %>% select(network, T_inter_blocks_mean, T_inter_aligned_mean,
                          T_inter_coherence), by = "network") %>%
  ## singles have no module pairs -> no inter-module coupling
  mutate(across(c(T_inter_blocks_mean, T_inter_aligned_mean),
                ~ifelse(n_modules == 1, 0, .x)))

INPUTS <- c(T_network            = "network T_s",
            T_network_intra_sum  = "sum of module T_s",
            T_inter_blocks_mean  = "inter-module T_s (team-aware, mean)",
            T_inter_aligned_mean = "inter-module alignment (signed)",
            T_inter_sum          = "inter-module T_s sum (old defn)",
            impurity             = "impurity",
            density              = "density",
            intra_density        = "intra-team density",
            inter_density        = "inter-team density",
            n_nodes              = "network size (nodes)")

## ---- generic correlation + plot -----------------------------------------
corr_tab <- function(dat, inputs, outputs, group = "n_modules") {
  bind_rows(lapply(sort(unique(dat[[group]])), function(k) {
    d <- dat[dat[[group]] == k, ]
    bind_rows(lapply(names(inputs), function(ix) bind_rows(lapply(names(outputs), function(oy) {
      x <- d[[ix]]; y <- d[[outputs[[oy]]]]
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0)
        return(data.frame(grp = k, input = ix, output = oy, r = NA_real_, p = NA_real_, n = sum(ok)))
      ct <- suppressWarnings(cor.test(x[ok], y[ok]))
      data.frame(grp = k, input = ix, output = oy, r = unname(ct$estimate), p = ct$p.value, n = sum(ok))
    }))))
  }))
}

draw <- function(tab, inputs, outputs, file, title, subtitle, caption, w, h) {
  nlab <- tab %>% group_by(grp) %>% summarise(n = max(n), .groups = "drop") %>%
    mutate(lab = sprintf("%s module%s  (n = %d)", grp, ifelse(grp == 1, "", "s"), n))
  t2 <- tab %>% left_join(nlab[, c("grp", "lab")], by = "grp") %>%
    mutate(lab = factor(lab, levels = nlab$lab),
           input = factor(inputs[input], levels = inputs),
           output = factor(output, levels = rev(names(outputs))),
           txt = ifelse(is.na(r), "", sprintf("%.2f%s", r, ifelse(!is.na(p) & p < 0.05, "*", ""))))
  p <- ggplot(t2, aes(input, output, fill = r)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    geom_text(aes(label = txt, colour = !is.na(r) & abs(r) > 0.6), size = 2.4, show.legend = FALSE) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), na.value = INK) +
    scale_fill_gradient2(low = BLUE, mid = GREY, high = RED, midpoint = 0, limits = c(-1, 1),
                         na.value = "grey94", name = "Pearson r", breaks = c(-1, -0.5, 0, 0.5, 1),
                         guide = guide_colourbar(barwidth = unit(0.32, "cm"),
                                                 barheight = unit(3.2, "cm"), ticks.colour = "white")) +
    facet_wrap(~lab, nrow = 1) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, caption = caption) +
    thm() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.72)),
                  axis.text.y = element_text(size = rel(0.78)),
                  legend.position = "right", legend.direction = "vertical")
  ggsave(file.path(outdir, file), p, width = w, height = h, dpi = 400, bg = "white")
  cat(sprintf("  %-32s %3d cells, %3d significant\n", file, sum(!is.na(t2$r)), sum(t2$p < 0.05, na.rm = TRUE)))
}

FOOT <- paste("Inputs at l_max = 10. Phenotypic scores and frustration: Ising steady states.",
              "PCA metrics: RACIPE continuous sampling. Inter-module T_s is the TEAM-AWARE block metric;",
              "the old whole-block definition is shown alongside for continuity and cancels team-aligned coupling to ~0.")
SMALL <- paste("* marks p < 0.05. Sample sizes are small at network grain: |r| > 0.81 is needed at n = 6, |r| > 0.44 at n = 20.",
               "At n = 1 the inter-module inputs are identically 0 (no module pairs), so those columns are blank.")

## ===== 1. PC and other variant-independent metrics ========================
PCOUT <- c("PC1 variance" = "PC1_var", "PCs to 90%" = "NumPC_90",
           "PCA entropy" = "PCA_entropy", "PCA entropy (norm)" = "PCA_entropy_norm",
           "frustration" = "frustration", "number of states" = "n_states")
draw(corr_tab(D, INPUTS, PCOUT), INPUTS, PCOUT, sprintf("corr_pc_lmax%d.png", LMAX),
     "PC and variant-independent landscape metrics",
     "outputs that do not depend on the phenotypic-score definition, so this matrix is shared by every variant",
     cap(FOOT, SMALL), 18.5, 4.6)
write.csv(corr_tab(D, INPUTS, PCOUT), file.path(outdir, sprintf("corr_pc_lmax%d.csv", LMAX)), row.names = FALSE)

## ===== 2-3. module-level P_s, mean AND variance ===========================
VLAB <- c(infl = "influence-weighted", state = "state-based")
for (v in names(VLAB)) {
  ps <- c("dominant", "top30", "top60")
  OUT <- setNames(c(sprintf("P_s_%s_%s__mean", v, ps), sprintf("P_s_%s_%s__var", v, ps)),
                  c(sprintf("MEAN over modules (%s)", ps), sprintf("VARIANCE over modules (%s)", ps)))
  tb <- corr_tab(D, INPUTS, OUT)
  draw(tb, INPUTS, OUT, sprintf("corr_module_%s_lmax%d.png", v, LMAX),
       sprintf("Module-level P_s -- %s", VLAB[[v]]),
       "rolled up to the network two ways: the MEAN over a network's modules (reproduces the published rows) and the VARIANCE (how unevenly the modules score)",
       cap(FOOT, "Variance is undefined at n = 1 (one module) and is computed over only 2 values at n = 2, so read those columns cautiously.", SMALL),
       18.5, 4.8)
  write.csv(tb, file.path(outdir, sprintf("corr_module_%s_lmax%d.csv", v, LMAX)), row.names = FALSE)
}

## ===== 4-5. network-level P_s =============================================
for (v in names(VLAB)) {
  ps <- c("dominant", "top30", "top60")
  OUT <- setNames(sprintf("P_s_net_%s_%s", v, ps), sprintf("network P_s (%s)", ps))
  tb <- corr_tab(D, INPUTS, OUT)
  draw(tb, INPUTS, OUT, sprintf("corr_network_%s_lmax%d.png", v, LMAX),
       sprintf("Network-level P_s -- %s", VLAB[[v]]),
       "the network's own two teams on the core influence matrix; no module roll-up is involved",
       cap(FOOT, SMALL), 18.5, 3.9)
  write.csv(tb, file.path(outdir, sprintf("corr_network_%s_lmax%d.csv", v, LMAX)), row.names = FALSE)
}

## ===== 6. TRUE module grain ==============================================
## One row per (network, module). Both sides are module-level, so there is no
## roll-up and no pseudo-replication.
io <- TA %>% group_by(network, module_from) %>%
        summarise(inter_out = sum(T_inter_blocks), inter_out_old = sum(T_inter_old), .groups = "drop") %>%
        rename(module = module_from) %>%
      full_join(TA %>% group_by(network, module_to) %>%
        summarise(inter_in = sum(T_inter_blocks), inter_in_old = sum(T_inter_old), .groups = "drop") %>%
        rename(module = module_to), by = c("network", "module")) %>%
      mutate(inter_total = inter_in + inter_out, inter_total_old = inter_in_old + inter_out_old)

msize <- read.csv("Metrics/module_properties.csv") %>% filter(l_max == LMAX) %>% select(module, mod_nodes = n_nodes)
G <- M %>% left_join(io, by = c("network", "module")) %>% left_join(msize, by = "module") %>%
  mutate(across(c(inter_in, inter_out, inter_total, inter_total_old),
                ~ifelse(n_modules == 1 | is.na(.x), 0, .x)))

GIN <- c(T_module_intra  = "this module's intra-team T_s",
         inter_in        = "inter-T_s INTO this module (team-aware)",
         inter_out       = "inter-T_s OUT of this module (team-aware)",
         inter_total     = "inter-T_s TOTAL for this module (team-aware)",
         inter_total_old = "inter-T_s total (old defn)",
         mod_nodes       = "module size (nodes)")
GOUT <- setNames(c(sprintf("P_s_infl_%s", c("dominant","top30","top60")),
                   sprintf("P_s_state_%s", c("dominant","top30","top60"))),
                 c(sprintf("module P_s, influence (%s)", c("dominant","top30","top60")),
                   sprintf("module P_s, state (%s)", c("dominant","top30","top60"))))
tb <- corr_tab(G, GIN, GOUT)
draw(tb, GIN, GOUT, sprintf("corr_modulegrain_lmax%d.png", LMAX),
     "Module grain: each module's own structure against its own P_s",
     "one row per (network, module) -- no roll-up, no repeated predictor, and 2-10x the sample size of the network-grain matrices",
     cap("Inputs at l_max = 10. Ising steady states.",
         "Inter-module coupling is summed over the module's ordered pairs using the TEAM-AWARE block metric; the old definition is shown for continuity.",
         "At n = 1 (the Comb_<M>_<M> self-pairs) there are no module pairs, so the inter-module columns are identically 0 and blank.",
         "* marks p < 0.05."), 17.0, 4.6)
write.csv(tb, file.path(outdir, sprintf("corr_modulegrain_lmax%d.csv", LMAX)), row.names = FALSE)

cat("\nnetwork-grain rows per module count:\n"); print(table(D$n_modules))
cat("module-grain rows per module count:\n");   print(table(G$n_modules))
cat("\ndone.\n")
