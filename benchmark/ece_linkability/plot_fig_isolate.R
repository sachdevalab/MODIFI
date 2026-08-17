#!/usr/bin/env Rscript
# Polished ECE linkage-ability figure -- ISOLATE dataset, plasmids and viruses POOLED.
# Reads the source-data CSVs from plot_fig_isolate.py. 1x2 layout:
#   a. ECE length distribution (density), linked vs unlinked
#   b. Host-motif density (recognition-site occ/kb): host chromosome vs linked vs unlinked ECE,
#      Welch t-test significance marked with stars.
# Run with:  /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_fig_isolate.R
suppressMessages({ library(ggplot2); library(ggsignif); library(patchwork) })

OUT <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_linkability"
A <- read.csv(file.path(OUT, "fig_isolate_linkability_a_length_sourcedata.csv"))
B <- read.csv(file.path(OUT, "fig_isolate_linkability_b_density_sourcedata.csv"))

C_HOST <- "#0072B2"; C_LINK <- "#009E73"; C_UNLINK <- "#D55E00"; GREY <- "#8A8A8A"
star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                    ifelse(p < 5e-2, "*", "ns"))))
base_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0),
        axis.text = element_text(color = "black"),
        legend.background = element_rect(fill = "#FFFFFFAA", color = NA),
        legend.key = element_blank(), plot.margin = margin(6, 10, 6, 6))

# ---------- a. length distribution ----------
A$group <- factor(A$group, c("linked", "unlinked"))
A$log10len <- log10(A$ece_len_bp)
med <- tapply(A$log10len, A$group, median)
n_lk <- sum(A$group == "linked"); n_ul <- sum(A$group == "unlinked")
pa <- t.test(log10len ~ group, data = A)$p.value
brks <- 3:6; lbls <- c(expression(10^3), expression(10^4), expression(10^5), expression(10^6))
pA <- ggplot(A, aes(log10len, fill = group, color = group)) +
  geom_density(alpha = 0.45, linewidth = 0.7, adjust = 1.1) +
  geom_vline(xintercept = med[["linked"]], color = C_LINK, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = med[["unlinked"]], color = C_UNLINK, linetype = "dashed", linewidth = 0.7) +
  scale_fill_manual(values = c(linked = C_LINK, unlinked = C_UNLINK),
                    labels = c(sprintf("linked (n=%d)", n_lk), sprintf("unlinked (n=%d)", n_ul)),
                    name = NULL) +
  scale_color_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), guide = "none") +
  scale_x_continuous(breaks = brks, labels = lbls) +
  guides(fill = guide_legend(override.aes = list(color = NA, alpha = 0.6))) +
  labs(x = "ECE length (bp)", y = "density",
       title = sprintf("a  Length distribution  (%s, p = %.1e)", star(pa), pa)) +
  base_theme +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.key.size = unit(0.85, "lines"), legend.text = element_text(size = 10))

# ---------- b. motif density ----------
lev <- c("host_chromosome", "linked_ECE", "unlinked_ECE")
B$group <- factor(B$group, lev)
ns <- table(B$group)
xlabs <- c(sprintf("host\nchromosome\n(n=%d)", ns[["host_chromosome"]]),
           sprintf("linked\nECE\n(n=%d)", ns[["linked_ECE"]]),
           sprintf("unlinked\nECE\n(n=%d)", ns[["unlinked_ECE"]]))
cols <- c(host_chromosome = C_HOST, linked_ECE = C_LINK, unlinked_ECE = C_UNLINK)
g <- function(k) B$density_per_kb[B$group == k]
p_hl <- t.test(g("host_chromosome"), g("linked_ECE"))$p.value
p_lu <- t.test(g("linked_ECE"), g("unlinked_ECE"))$p.value
p_hu <- t.test(g("host_chromosome"), g("unlinked_ECE"))$p.value
top <- as.numeric(max(B$density_per_kb, na.rm = TRUE))   # above all jittered points/whiskers
y_lo <- top * 1.08; y_hi <- top * 1.26
pB <- ggplot(B, aes(group, density_per_kb, fill = group)) +
  geom_jitter(aes(color = group), width = 0.11, size = 0.5, alpha = 0.22) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.85, linewidth = 0.5, color = "grey20") +
  geom_signif(annotations = c(star(p_hl), star(p_lu), star(p_hu)),
              xmin = c(1, 2, 1), xmax = c(2, 3, 3),
              y_position = c(y_lo, y_lo, y_hi),
              tip_length = 0.008, textsize = 5, vjust = 0.2) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_color_manual(values = cols, guide = "none") +
  scale_x_discrete(labels = xlabs) +
  coord_cartesian(ylim = c(0, top * 1.42)) +
  labs(x = NULL, y = "motif density (recognition-site occ / kb)",
       title = "b  Host-motif density") +
  base_theme

fig <- pA + pB + plot_layout(widths = c(1, 1)) +
  plot_annotation(title = "Isolate ECEs (1,120 evaluable plasmids and viruses)",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5)))

ggsave(file.path(OUT, "fig_isolate_linkability.pdf"), fig, width = 10.5, height = 4.6, device = cairo_pdf)
ggsave(file.path(OUT, "fig_isolate_linkability.png"), fig, width = 10.5, height = 4.6, dpi = 300)
cat(sprintf("wrote fig_isolate_linkability.pdf/.png (pooled; length %s p=%.2e)\n", star(pa), pa))
