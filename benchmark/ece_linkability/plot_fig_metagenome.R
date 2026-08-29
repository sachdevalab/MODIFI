#!/usr/bin/env Rscript
# Combined 3-panel metagenome figure (strict ECE set, 64 run2 metagenomes).
#   a  length distribution of linked vs unlinked strict ECEs (pooled plasmid+virus)
#   b  paired scatter: linked host motif-set density (x) vs same set on the ECE (y), per linked ECE
#   c  sample-level metagenome-motif density: linked host vs linked ECE vs unlinked ECE
#      (one point per sample per group; paired Wilcoxon across samples, controls sample complexity)
# Run with: /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_fig_metagenome.R
suppressMessages({ library(ggplot2); library(ggsignif); library(patchwork) })

OUT <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_linkability"
cargs <- commandArgs(trailingOnly = TRUE)
PRE <- if (length(cargs) >= 1) cargs[1] else ""   # "loose_" for the loose-linkage version
A <- read.csv(file.path(OUT, sprintf("fig_metagenome_%sa_length_sourcedata.csv", PRE)))
B <- read.csv(file.path(OUT, sprintf("fig_metagenome_%sb_scatter_sourcedata.csv", PRE)))
Cc <- read.csv(file.path(OUT, sprintf("fig_metagenome_%sc_density_sourcedata.csv", PRE)))

C_HOST <- "#0072B2"; C_LINK <- "#009E73"; C_UNLINK <- "#D55E00"; GREY <- "#8A8A8A"
star <- function(p) ifelse(is.na(p), "ns", ifelse(p < 1e-3, "***", ifelse(p < 1e-2, "**",
                    ifelse(p < 5e-2, "*", "ns"))))
base_theme <- theme_classic(base_size = 12) +
  theme(plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.tag = element_text(face = "bold", size = 15),
        axis.text = element_text(color = "black"),
        legend.background = element_rect(fill = "#FFFFFFAA", color = NA),
        legend.key = element_blank(), plot.margin = margin(6, 10, 6, 6))

# ---------- a. length distribution linked vs unlinked ----------
A$grp <- factor(ifelse(A$linked == "True" | A$linked == TRUE, "linked", "unlinked"),
                c("linked", "unlinked"))
A$log10len <- log10(A$mge_len)
med <- tapply(A$log10len, A$grp, median)
n_lk <- sum(A$grp == "linked"); n_ul <- sum(A$grp == "unlinked")
pa <- t.test(log10len ~ grp, data = A)$p.value
brks <- 3:6; lbls <- c(expression(10^3), expression(10^4), expression(10^5), expression(10^6))
pA <- ggplot(A, aes(log10len, fill = grp, color = grp)) +
  geom_density(alpha = 0.45, linewidth = 0.7, adjust = 1.1) +
  geom_vline(xintercept = med[["linked"]], color = C_LINK, linetype = "dashed", linewidth = 0.7) +
  geom_vline(xintercept = med[["unlinked"]], color = C_UNLINK, linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = 3.62, y = 1.5, hjust = 0, size = 3.3,
           label = sprintf("%s  p = %.1e", star(pa), pa)) +
  scale_fill_manual(values = c(linked = C_LINK, unlinked = C_UNLINK),
                    labels = c("linked", "unlinked"),
                    name = NULL) +
  scale_color_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), guide = "none") +
  scale_x_continuous(breaks = brks, labels = lbls) +
  guides(fill = guide_legend(override.aes = list(color = NA, alpha = 0.6))) +
  labs(x = "ECE length (bp)", y = "density",
       title = sprintf("b  Length distribution  (%s, p = %.1e)", star(pa), pa)) +
  base_theme +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.key.size = unit(0.85, "lines"), legend.text = element_text(size = 9))

# ---------- b. paired scatter: host motif density vs ECE ----------
B$type <- factor(B$type, c("plasmid", "virus"))
cols <- c(plasmid = C_LINK, virus = C_UNLINK)
lim <- as.numeric(quantile(c(B$host_density, B$ece_density), 0.975, na.rm = TRUE))  # zoom to the bulk
below <- function(t) { s <- B[B$type == t, ]; mean(s$ece_density < s$host_density) * 100 }
pw <- function(t) { s <- B[B$type == t, ]; wilcox.test(s$host_density, s$ece_density, paired = TRUE)$p.value }
n_pl <- sum(B$type == "plasmid"); n_vi <- sum(B$type == "virus")
labB <- sprintf(paste0("paired Wilcoxon (ECE vs host)\n",
                       "plasmid: %.0f%% below, p = %.1e\n",
                       "virus: %.0f%% below, p = %.1e"),
                below("plasmid"), pw("plasmid"), below("virus"), pw("virus"))
pB <- ggplot(B, aes(host_density, ece_density, color = type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_point(size = 1.1, alpha = 0.55) +
  annotate("text", x = lim * 0.14, y = lim, hjust = 0, vjust = 1, size = 2.9, label = labB) +
  scale_color_manual(values = cols, name = NULL,
                     labels = c("plasmid", "virus")) +
  coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
  labs(x = "linked host motif density (occ / kb)", y = "ECE motif density (occ / kb)",
       title = "a  Host motif density on linked ECEs") +
  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
  base_theme + theme(legend.position = c(0.98, 0.03), legend.justification = c(1, 0))

# ---------- d & e. modification site number / frequency, linked vs unlinked ----------
luplot <- function(col, ylab, title, logy) {
  d <- A[is.finite(A[[col]]), ]
  p <- wilcox.test(d[[col]] ~ d$grp)$p.value
  n_lk <- sum(d$grp == "linked"); n_ul <- sum(d$grp == "unlinked")
  top <- as.numeric(quantile(d[[col]], 0.985, na.rm = TRUE))
  segy <- if (logy) 10^(log10(top) + 0.15) else top * 1.05
  txty <- if (logy) 10^(log10(top) + 0.32) else top * 1.13
  ylimtop <- if (logy) 10^(log10(top) + 0.55) else top * 1.25
  g <- ggplot(d, aes(grp, .data[[col]], fill = grp)) +
    geom_jitter(aes(color = grp), width = 0.12, size = 0.4, alpha = 0.18) +
    geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.85, linewidth = 0.5, color = "grey20") +
    annotate("segment", x = 1, xend = 2, y = segy, yend = segy, linewidth = 0.4) +
    annotate("text", x = 1.5, y = txty, label = sprintf("%s  p = %.1e", star(p), p), size = 3.6) +
    scale_fill_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), guide = "none") +
    scale_color_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), guide = "none") +
    scale_x_discrete(labels = c(sprintf("linked\n(%d)", n_lk), sprintf("unlinked\n(%d)", n_ul))) +
    labs(x = NULL, y = ylab, title = title) +
    base_theme
  if (logy) g <- g + scale_y_log10()
  g + coord_cartesian(ylim = c(if (logy) NA else 0, ylimtop))
}
p_modn <- luplot("n_mod_sites", "modification sites (score >= 30)", "e  Modification site number", TRUE)
p_modf <- luplot("mod_density_per_kb", "modified sites / kb", "d  Modified site density", FALSE)

# ---------- e. modification frequency by length bin (length-controlled) ----------
dc <- A[is.finite(A$mod_density_per_kb) & is.finite(A$mge_len), ]
dc$lk <- factor(ifelse(dc$linked == TRUE | dc$linked == "True", "linked", "unlinked"),
                c("linked", "unlinked"))
edges_l <- unique(quantile(dc$mge_len, seq(0, 1, 0.2)))
dc$bin <- cut(dc$mge_len, breaks = edges_l, include.lowest = TRUE)
lv <- levels(dc$bin)
binlab <- sapply(strsplit(gsub("[][()]", "", lv), ","), function(x)
  sprintf("%g-%g\nkb", round(as.numeric(x[1]) / 1000), round(as.numeric(x[2]) / 1000)))
mlm <- lm(mod_density_per_kb ~ lk + log10(mge_len), data = dc)
p_link <- summary(mlm)$coefficients["lkunlinked", 4]
binp <- sapply(lv, function(b) {
  g <- dc[dc$bin == b, ]; L <- g$mod_density_per_kb[g$lk == "linked"]; U <- g$mod_density_per_kb[g$lk == "unlinked"]
  if (length(L) >= 3 && length(U) >= 3) wilcox.test(L, U)$p.value else NA
})
ptxt <- ifelse(is.na(binp) | binp >= 0.05, "n.s.", sprintf("p=%.3f", binp))
ytop <- as.numeric(quantile(dc$mod_density_per_kb, 0.99, na.rm = TRUE))
p_lenctrl <- ggplot(dc, aes(bin, mod_density_per_kb, fill = lk)) +
  geom_boxplot(outlier.shape = NA, width = 0.7, linewidth = 0.3, position = position_dodge(0.8)) +
  annotate("text", x = seq_along(lv), y = ytop * 1.03, label = ptxt, size = 3) +
  scale_fill_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), name = NULL) +
  scale_x_discrete(labels = binlab) +
  coord_cartesian(ylim = c(0, ytop * 1.1)) +
  labs(x = "ECE length bin", y = "modified sites / kb",
       title = "f  Modified site density by length bin",
       subtitle = sprintf("length-adjusted (OLS): %s (p = %.2f)",
                          ifelse(p_link < 0.05, "sig.", "n.s."), p_link)) +
  base_theme + theme(legend.position = "top", legend.justification = "left",
                     plot.subtitle = element_blank())

# ---------- e & f. any metagenome motif: linked vs unlinked ECE, per sample ----------
mk_scatter <- function(valcol, axlab, title) {
  med <- aggregate(reformulate(c("sample", "group"), valcol), data = Cc, FUN = median)
  ww <- as.data.frame(tapply(med[[valcol]], list(med$sample, med$group), function(x) x[1]))
  d <- ww[complete.cases(ww$linked_ECE, ww$unlinked_ECE), ]
  p <- wilcox.test(d$linked_ECE, d$unlinked_ECE, paired = TRUE)$p.value
  gt <- mean(d$linked_ECE > d$unlinked_ECE) * 100
  lab <- sprintf("linked > unlinked: %.0f%%\npaired %s, p = %.1e", gt, star(p), p)
  lo <- min(d$linked_ECE, d$unlinked_ECE, na.rm = TRUE); hi <- max(d$linked_ECE, d$unlinked_ECE, na.rm = TRUE)
  ggplot(d, aes(linked_ECE, unlinked_ECE)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_point(size = 1.5, alpha = 0.75, color = "#4477AA") +
    annotate("text", x = lo, y = hi, hjust = 0, vjust = 1, size = 2.7, label = lab) +
    scale_x_log10() + scale_y_log10() + coord_equal(xlim = c(lo, hi), ylim = c(lo, hi)) +
    labs(x = sprintf("linked ECE %s", axlab), y = sprintf("unlinked ECE %s", axlab), title = title) +
    base_theme
}
p_amd <- mk_scatter("density_per_kb", "density (occ/kb)", "c  Metagenome-motif density")
p_amo <- mk_scatter("n_occ", "occurrences (count)", "d  Metagenome-motif occurrences")

# ---------- f. metagenome-motif density: WITHIN-SAMPLE, by length bin ----------
# Length-binned AND same-sample-controlled: within each length bin, for every sample that has BOTH
# linked and unlinked ECEs in that bin, take the within-sample log2 ratio of their median densities.
# Boxes at 0 => equal; only samples with both groups in a bin contribute ("otherwise not compared").
cf <- Cc[Cc$group %in% c("linked_ECE", "unlinked_ECE") & is.finite(Cc$density_per_kb) &
         is.finite(Cc$ece_len) & Cc$density_per_kb > 0, ]
cf$lk <- ifelse(cf$group == "linked_ECE", "linked", "unlinked")
edges_f <- unique(quantile(cf$ece_len, seq(0, 1, 0.2)))
cf$bin <- cut(cf$ece_len, breaks = edges_f, include.lowest = TRUE)
lvf <- levels(cf$bin)
agg <- aggregate(density_per_kb ~ sample + bin + lk, data = cf, FUN = median)
rlist <- list()
for (b in lvf) {
  gb <- agg[agg$bin == b, ]
  L <- gb$density_per_kb[gb$lk == "linked"]; names(L) <- gb$sample[gb$lk == "linked"]
  U <- gb$density_per_kb[gb$lk == "unlinked"]; names(U) <- gb$sample[gb$lk == "unlinked"]
  cm <- intersect(names(L), names(U))
  if (length(cm) >= 3) rlist[[b]] <- data.frame(bin = b, lr = log2(L[cm] / U[cm]))
}
ratios <- do.call(rbind, rlist)
ratios$bin <- factor(ratios$bin, lvf)
binlabf <- sapply(strsplit(gsub("[][()]", "", lvf), ","), function(x)
  sprintf("%g-%g\nkb", round(as.numeric(x[1]) / 1000), round(as.numeric(x[2]) / 1000)))
binpf <- sapply(lvf, function(b) { x <- ratios$lr[ratios$bin == b]; if (length(x) >= 3) wilcox.test(x)$p.value else NA })
ptxtf <- ifelse(is.na(binpf) | binpf >= 0.05, "n.s.", sprintf("p=%.3f", binpf))
# paired long form: for each length bin, the samples that have BOTH groups (connected by a line)
paired_long <- do.call(rbind, lapply(lvf, function(b) {
  gb <- agg[agg$bin == b, ]; L <- gb[gb$lk == "linked", ]; U <- gb[gb$lk == "unlinked", ]
  cm <- intersect(L$sample, U$sample)
  if (length(cm) >= 3) rbind(
    data.frame(bin = b, sample = cm, lk = "linked", y = L$density_per_kb[match(cm, L$sample)]),
    data.frame(bin = b, sample = cm, lk = "unlinked", y = U$density_per_kb[match(cm, U$sample)])) else NULL
}))
paired_long$bin <- factor(paired_long$bin, lvf)
paired_long$lk <- factor(paired_long$lk, c("linked", "unlinked"))
strip <- setNames(sapply(strsplit(gsub("[][()]", "", lvf), ","), function(x)
  sprintf("%g-%gkb", round(as.numeric(x[1]) / 1000), round(as.numeric(x[2]) / 1000))), lvf)
statf <- data.frame(bin = factor(lvf, lvf), lab = ptxtf, y = max(paired_long$y, na.rm = TRUE))
p_lenctrl2 <- ggplot(paired_long, aes(lk, y)) +
  geom_line(aes(group = sample), color = "grey70", alpha = 0.5, linewidth = 0.25) +
  geom_point(aes(color = lk), size = 1, alpha = 0.85) +
  geom_text(data = statf, aes(x = 1.5, y = y, label = lab), size = 2.7, vjust = 1, inherit.aes = FALSE) +
  facet_wrap(~bin, nrow = 1, labeller = labeller(bin = strip)) +
  scale_color_manual(values = c(linked = C_LINK, unlinked = C_UNLINK), guide = "none") +
  scale_y_log10() + scale_x_discrete(labels = c("linked", "unlinked")) +
  labs(x = NULL, y = "metagenome-motif density (occ/kb)",
       title = "e  Metagenome-motif density, within-sample x length",
       subtitle = "lines = same sample (both groups in bin); signed-rank per bin") +
  base_theme + theme(plot.subtitle = element_blank(),
                     axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
                     panel.spacing = unit(2, "pt"))

mode_lab <- if (PRE == "loose_") "loose linkage: score>0.5 & spec<0.01" else "high-confidence linkage (curated)"
fig <- pB + pA + p_amd + p_modf + p_lenctrl2 + p_lenctrl +
  plot_layout(ncol = 3, widths = c(1, 2, 1.2), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 15)))
ggsave(file.path(OUT, sprintf("fig_metagenome_%slinkability.pdf", PRE)), fig, width = 20, height = 9, device = cairo_pdf)
ggsave(file.path(OUT, sprintf("fig_metagenome_%slinkability.png", PRE)), fig, width = 20, height = 9, dpi = 190)
cat(sprintf("wrote fig_metagenome_%slinkability.* (%d linked; a length %s; b plasmid %s virus %s)\n",
            PRE, n_lk, star(pa), star(pw("plasmid")), star(pw("virus"))))
