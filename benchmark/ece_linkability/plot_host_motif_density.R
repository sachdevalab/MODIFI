#!/usr/bin/env Rscript
# Standalone paired scatter (no title): host chromosome motif density (x) vs ECE motif density (y),
# strict isolate ECE set (known ECE-host linkage). Points below the y=x line are ECEs depleted of
# host recognition motifs relative to their own host chromosome (restriction-avoidance). Paired
# Wilcoxon signed-rank p-values annotated for plasmids and viruses.
# Run with: /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_host_motif_density.R
suppressMessages({ library(ggplot2) })

OUT <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_linkability"
STEM <- "fig_host_motif_density"
B <- read.csv(file.path(OUT, paste0(STEM, "_sourcedata.csv")))
B$type <- factor(B$type, c("plasmid", "virus"))
cols <- c(plasmid = "#009E73", virus = "#D55E00")

lim <- max(B$host_density, B$ece_density, na.rm = TRUE) * 1.02
pw <- function(t) {
  s <- B[B$type == t, ]
  wilcox.test(s$host_density, s$ece_density, paired = TRUE)$p.value
}
below <- function(t) { s <- B[B$type == t, ]; mean(s$ece_density < s$host_density) * 100 }
p_pl <- pw("plasmid"); p_vi <- pw("virus")
n_pl <- sum(B$type == "plasmid"); n_vi <- sum(B$type == "virus")
lab <- sprintf(paste0("plasmid (n=%d): %.0f%% below y=x, Wilcoxon paired p = %.1e\n",
                      "virus (n=%d): %.0f%% below y=x, Wilcoxon paired p = %.1e"),
               n_pl, below("plasmid"), p_pl, n_vi, below("virus"), p_vi)

p <- ggplot(B, aes(host_density, ece_density, color = type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_point(size = 1.1, alpha = 0.5) +
  annotate("text", x = lim * 0.05, y = lim, hjust = 0, vjust = 1, size = 3.6, label = lab) +
  scale_color_manual(values = cols, name = NULL,
                     labels = c(sprintf("plasmid (n=%d)", n_pl), sprintf("virus (n=%d)", n_vi))) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
  labs(x = "host chromosome motif density (recognition-site occ / kb)",
       y = "ECE motif density (recognition-site occ / kb)") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.98, 0.03), legend.justification = c(1, 0),
        plot.margin = margin(8, 12, 8, 8))

ggsave(file.path(OUT, paste0(STEM, ".pdf")), p, width = 5.8, height = 5.6, device = cairo_pdf)
ggsave(file.path(OUT, paste0(STEM, ".png")), p, width = 5.8, height = 5.6, dpi = 300)
cat(sprintf("wrote %s.pdf/.png  (plasmid p=%.2e, virus p=%.2e)\n", STEM, p_pl, p_vi))
