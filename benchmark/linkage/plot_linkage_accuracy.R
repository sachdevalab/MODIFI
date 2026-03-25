library(ggplot2)
library(dplyr)
library(patchwork)

fig_dir <- "../../tmp/figures/link_accuracy/"

curve_df  <- read.csv(file.path(fig_dir, "linkage_curve.csv"))
cutoff_df <- read.csv(file.path(fig_dir, "linkage_cutoffs.csv"))

# ── Explicit label positions (hand-tuned per panel) ────────────────────────────
# recall panel: y range 0–1.08
cutoff_df <- cutoff_df %>%
  mutate(
    lx_recall = c( 0.15,  1.12,  1.45),
    ly_recall = c( 0.34,  0.74,  0.96),
    lx_prec   = c( 0.15,  0.15,  1.45),
    ly_prec   = c( 0.975, 0.945, 0.906)
  )

# ── Theme ──────────────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "grey40"),
    axis.ticks       = element_line(colour = "grey40"),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0),
    axis.title       = element_text(size = 11),
    legend.position  = "none"
  )

curve_col <- "#2166AC"
point_col <- "#C0392B"
total_fp  <- max(curve_df$fp)

# ── Panel builder ──────────────────────────────────────────────────────────────
make_panel <- function(y_col, lx_col, ly_col, ylab, ylim, ybreaks) {
  ggplot() +
    geom_step(
      data = curve_df,
      aes(x = fp, y = .data[[y_col]]),
      colour = curve_col, linewidth = 1.2, direction = "hv"
    ) +
    geom_vline(
      xintercept = c(1, 2),
      linetype = "dashed", colour = "grey70", linewidth = 0.35
    ) +
    # connector segments: point → label
    geom_segment(
      data = cutoff_df,
      aes(x = fp, y = .data[[y_col]],
          xend = .data[[lx_col]], yend = .data[[ly_col]]),
      colour = "grey55", linewidth = 0.35
    ) +
    # label boxes
    geom_label(
      data = cutoff_df,
      aes(x = .data[[lx_col]], y = .data[[ly_col]], label = label),
      hjust = 0, size = 3.1, colour = "grey15", fill = "white",
      label.padding = unit(0.2, "lines"), label.size = 0.2
    ) +
    # points on top
    geom_point(
      data = cutoff_df,
      aes(x = fp, y = .data[[y_col]]),
      shape = 21, size = 4, fill = "white",
      colour = point_col, stroke = 1.5
    ) +
    scale_x_continuous(
      breaks = 0:total_fp,
      limits = c(-0.1, total_fp + 0.1),
      expand = expansion(add = 0.05)
    ) +
    scale_y_continuous(
      limits = ylim,
      breaks = ybreaks,
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(add = 0)
    ) +
    labs(x = "Number of false positives", y = ylab) +
    base_theme
}

p_recall <- make_panel(
  y_col   = "recall",
  lx_col  = "lx_recall", ly_col = "ly_recall",
  ylab    = "Recall",
  ylim    = c(0, 1.10),
  ybreaks = seq(0, 1, 0.25)
)

p_precision <- make_panel(
  y_col   = "precision",
  lx_col  = "lx_prec", ly_col = "ly_prec",
  ylab    = "Precision",
  ylim    = c(0.88, 1.02),
  ybreaks = seq(0.88, 1.0, 0.04)
)

# ── Save ───────────────────────────────────────────────────────────────────────
p_combined <- p_recall | p_precision

ggsave(file.path(fig_dir, "linkage_accuracy.pdf"),
       p_combined, width = 10, height = 4.5)
ggsave(file.path(fig_dir, "linkage_recall.pdf"),
       p_recall, width = 5, height = 4.5)
ggsave(file.path(fig_dir, "linkage_precision.pdf"),
       p_precision, width = 5, height = 4.5)

cat("Saved to", fig_dir, "\n")
