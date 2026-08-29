#!/usr/bin/env Rscript
# Separate figure: the WHOLE metagenome motif repertoire ("any metagenome motif") on linked vs
# unlinked ECEs, per sample (both axes are ECEs in the same sample, so this controls sample
# complexity). Each point = one metagenome sample; x = per-sample median on linked ECEs, y = on
# unlinked ECEs. Two subpanels:
#   a  density (occurrences / kb)
#   b  absolute occurrences (count)
# The dashed line is y = x (equal on linked and unlinked ECEs).
# Run with: /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_fig_metagenome_anymotif.R
suppressMessages({ library(ggplot2); library(patchwork) })

OUT <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_linkability"
STEM <- "fig_metagenome_anymotif"
Cc <- read.csv(file.path(OUT, "fig_metagenome_c_density_sourcedata.csv"))
PT <- "#4477AA"
star <- function(p) ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**", ifelse(p < 5e-2, "*", "ns")))

mk_scatter <- function(valcol, axlab, title) {
  med <- aggregate(reformulate(c("sample", "group"), valcol), data = Cc, FUN = median)
  w <- as.data.frame(tapply(med[[valcol]], list(med$sample, med$group), function(x) x[1]))
  d <- w[complete.cases(w$linked_ECE, w$unlinked_ECE), ]
  p <- wilcox.test(d$linked_ECE, d$unlinked_ECE, paired = TRUE)$p.value
  above <- mean(d$unlinked_ECE > d$linked_ECE) * 100   # unlinked higher than linked
  lab <- sprintf("n = %d samples\nunlinked > linked in %.0f%%\npaired %s, p = %.1e",
                 nrow(d), above, star(p), p)
  lo <- min(d$linked_ECE, d$unlinked_ECE, na.rm = TRUE)
  hi <- max(d$linked_ECE, d$unlinked_ECE, na.rm = TRUE)
  ggplot(d, aes(linked_ECE, unlinked_ECE)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_point(size = 1.9, alpha = 0.75, color = PT) +
    annotate("text", x = lo, y = hi, hjust = 0, vjust = 1, size = 3.2, label = lab) +
    scale_x_log10() + scale_y_log10() +
    coord_equal(xlim = c(lo, hi), ylim = c(lo, hi)) +
    labs(x = sprintf("linked ECE %s", axlab), y = sprintf("unlinked ECE %s", axlab), title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
          axis.text = element_text(color = "black"))
}

pA <- mk_scatter("density_per_kb", "density (occ / kb, per sample)", "a  Density")
pB <- mk_scatter("n_occ", "occurrences (count, per sample)", "b  Absolute occurrences")
fig <- pA + pB +
  plot_annotation(title = "Any metagenome motif: linked vs unlinked ECE, per sample",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5)))
ggsave(file.path(OUT, paste0(STEM, ".pdf")), fig, width = 11, height = 5.6, device = cairo_pdf)
ggsave(file.path(OUT, paste0(STEM, ".png")), fig, width = 11, height = 5.6, dpi = 220)
cat(sprintf("wrote %s.* (linked vs unlinked ECE; density + absolute occ)\n", STEM))
