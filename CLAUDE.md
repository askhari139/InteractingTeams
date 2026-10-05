# Project notes for Claude

## Language preference

**Write analysis code in R.** This is the preferred language for this project —
use it for new scripts, metric calculations, and plotting unless I explicitly
ask for something else.

Useful context:
- The `funcsKishore` package provides the core primitives: `InfluenceMatrix()`,
  `getGsVec()`, `TopoToIntMat()`, `ComputePowerMatrix()`, `findClusterTeams()`,
  `peripherals()`, plus `theme_Publication()` / `theme_Heatmap()` for plots.
- `script.r` at the repo root is the existing scratch entry point.
- Note `InfluenceMatrix()` caches to `Influence/<net>_reducedInfl.csv` and
  `setwd()`s during the call, so run it with the working directory set to the
  folder holding the `.topo` file; pass `force = TRUE` to bypass the cache.

## Figure preparation

Follow the conventions of the regulatory-noise paper code:
`~/Desktop/PostDoc/Abhay_Lakshmi/Regulatory-noise-model/codes/R/`
(`figure_common.r`, `run_all_figures.r`, `Figure1.r`–`Figure6.r`). Check those files
before building a new figure type.

- **Libraries**: `funcsKishore`, `tidyverse`, `cowplot` (assembly), `patchwork`,
  `ggpubr` (`stat_cor`), `ggsignif` (significance brackets).
- **Script layout**: header comment with `PURPOSE` and a per-panel list; a
  `# ── CONFIG ──` block; shared helpers in a `figure_common.r` sourced by each
  `FigureN.r`. A panel function builds the plot, saves it with `ggsave()` if
  `outFile` is given, and returns `invisible(p)`.
- **Theme**: `theme_Publication()` on every panel (`theme_Heatmap()` for heatmaps).
  Common tweaks: `axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)`,
  `strip.text = element_text(size = rel(1.2))`, `legend.position = "top"` if the
  default bottom legend crowds the panel.
- **Output**: `.jpg` via `ggsave()` with sizes in inches (default dpi 300).
  Individual panels go to `figures/individual/` (about 5x4, 6x4, 5x5 or 6x6 in).
  Combined figures go to `figures/final/FigN.jpg` / `FigNS1.jpg`.
- **Assembly**: nested `cowplot::plot_grid()` with uppercase `labels = c("A", ...)`,
  `label_size = LS` (`LS = 20`), and `rel_widths` / `rel_heights` for layout.
  Final size is `FIG_WIDTH` x `FIG_HEIGHT` (typically 12–20 in x 8–20 in), falling
  back to script defaults via `if (exists("FIG_WIDTH")) FIG_WIDTH else <default>`.
  A master `run_all_figures.r` sets `LS`, `FIG_WIDTH` and `FIG_HEIGHT` per script
  and sources it.
- **Colour**: discrete categories use named `scales::hue_pal()` vectors, kept
  consistent across panels. Ordered or continuous variables use viridis
  (`scale_color_viridis_d/_c`; `trans = "log10"` for densities). Correlations use
  `scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
  midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1))`, with white tile borders
  (`linewidth = 0.5`) and `geom_text(size = 3)`.
- **Statistics**: mean ± SD via `stat_summary(fun.data = mean_sd, geom = "errorbar",
  width = 0.3)` plus mean points (`size = 2`), `position_dodge(0.6)`. Paired Wilcoxon
  results appear as stars through `geom_signif(manual = TRUE, tip_length = 0.01,
  size = 0.3, textsize = 2.6)`. Spearman uses `ggpubr::stat_cor(output.type = "text")`
  showing `R = x` and `" *"` when p < 0.05, not the exact p-value.
- **Marks**: scatter `alpha = 0.6, size = 1`; y = x reference line
  `geom_abline(color = "black", linewidth = 0.6)`.
