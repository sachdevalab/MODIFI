#!/usr/bin/env Rscript
# Combined per-site IPD-ratio distributions for the two RS II 5mC test cases (WGA control):
# O. crateris (6mA/4mC/5mC) + T. vulgaris (two 5mC motifs). Panels a-e, no in-plot annotations
# (all described in the caption). Input from Python-binned densities.
suppressMessages({library(ggplot2); library(dplyr)})

d   <- read.csv("/home/shuaiw/borg/revision/ipd_combined_bins.csv", stringsAsFactors = FALSE)
med <- read.csv("/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification/ipd_distributions_combined_sourcedata.csv",
                stringsAsFactors = FALSE)
out <- "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"

# per-panel colour for the motif curve (control is grey in every panel)
mcol <- c(a = "#C0392B", b = "#2E86C1", c = "#E59400", d = "#8E44AD", e = "#16A085")
d$panel <- factor(d$panel, levels = letters[1:5])

ctrl  <- d %>% filter(group == "control")
motif <- d %>% filter(group == "motif")
mmed  <- med %>% transmute(panel = factor(panel, levels = letters[1:5]), x = motif_median)

dummy <- data.frame(ipd_ratio = c(-1, -1), density = c(0, 0), series = "motif")   # off-canvas
p <- ggplot() +
  # grey control histogram (fill legend), motif as a coloured step outline (linetype legend)
  geom_col(data = ctrl,  aes(ipd_ratio, density, fill = "control"), width = 5/60, alpha = 0.5) +
  geom_step(data = motif, aes(ipd_ratio, density, colour = panel), linewidth = 0.8,
            direction = "mid") +
  geom_line(data = dummy, aes(ipd_ratio, density, linetype = series)) +
  geom_vline(xintercept = 1.0, linetype = "dotted", linewidth = 0.4) +
  geom_vline(data = mmed, aes(xintercept = x, colour = panel), linetype = "dashed", linewidth = 0.5) +
  geom_text(data = data.frame(panel = factor(letters[1:5], levels = letters[1:5])),
            aes(x = -Inf, y = Inf, label = panel), fontface = "bold", size = 5.5,
            hjust = -0.4, vjust = 1.3) +
  # strain (italic) + motif string, top-right of each panel
  geom_text(data = data.frame(
              panel = factor(letters[1:5], levels = letters[1:5]),
              lab = c('atop(italic("O. crateris"), "ATTAAT")',
                      'atop(italic("O. crateris"), "GGATCC")',
                      'atop(italic("O. crateris"), "GCNGC")',
                      'atop(italic("T. vulgaris"), "GGCC")',
                      'atop(italic("T. vulgaris"), "CCWGG")')),
            aes(x = Inf, y = Inf, label = lab), parse = TRUE, size = 3.6,
            hjust = 1.05, vjust = 1.15, lineheight = 0.9) +
  scale_colour_manual(values = mcol, guide = "none") +
  scale_fill_manual(name = NULL, values = c(control = "grey60")) +
  scale_linetype_manual(name = NULL, values = c(motif = "solid")) +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  scale_x_continuous("IPD ratio", limits = c(0, 5), expand = c(0.01, 0)) +
  ylab("density") +
  guides(fill     = guide_legend(order = 1, override.aes = list(alpha = 0.5)),
         linetype = guide_legend(order = 2, override.aes = list(colour = "grey20"))) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
        legend.position = c(0.86, 0.28),
        legend.spacing.y = unit(0, "pt"),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))

ggsave(file.path(out, "ipd_distributions_combined.pdf"), p, width = 9.6, height = 5.6)
ggsave(file.path(out, "ipd_distributions_combined.png"), p, width = 9.6, height = 5.6, dpi = 200)
message("saved ipd_distributions_combined.pdf/.png")
