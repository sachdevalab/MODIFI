#!/usr/bin/env Rscript
# plot_benchmark_robustness.R -- base-R integrated comprehensive benchmark figure.
# Panels: a complexity, b sequencing depth, c end-to-end de-novo (all recall/precision),
# d orphan-ECE false-positive control (ROC). Reads the Source Data CSVs so numbers match
# the Python figures exactly. Base graphics only (no ggplot). Panel titles are bare letters;
# all descriptive text lives in the caption.

BASE <- "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta"
SRC  <- file.path(BASE, "benchmark_summary", "benchmark_robustness_sourcedata.csv")
ROC  <- file.path(BASE, "orphan_fpr", "orphan_roc_pr_sourcedata.csv")
OUT  <- file.path(BASE, "benchmark_summary", "benchmark_robustness.pdf")
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREY <- "#EEEEEE"

d <- read.csv(SRC, stringsAsFactors = FALSE)
roc <- read.csv(ROC, stringsAsFactors = FALSE)
roc <- roc[order(roc$fpr_grid), ]
# AUROC of the mean ROC curve (trapezoid)
auroc <- sum(diff(roc$fpr_grid) * (head(roc$mean_tpr, -1) + tail(roc$mean_tpr, -1)) / 2)

metric_panel <- function(sub, xlab, letter, logx = FALSE, shade = NULL) {
  sub <- sub[order(sub$x), ]; x <- sub$x; xr <- range(x)
  plot(NA, xlim = if (logx) xr else xr + c(-0.04, 0.04) * diff(xr),
       ylim = c(0, 1.08), log = if (logx) "x" else "",
       xlab = xlab, ylab = "metric", main = "", xaxt = "n", las = 1)
  title(main = letter, adj = 0, font.main = 2, cex.main = 1.7, line = 0.6)
  axis(1, at = x, labels = x)
  grid(nx = NA, ny = NULL, col = "grey85", lty = 1)
  if (!is.null(shade)) rect(shade[1], -1, shade[2], 2, col = GREY, border = NA)
  box()
  eb <- function(y, ci, col, pch) {
    ok <- ci > 0
    if (any(ok)) arrows(x[ok], (y - ci)[ok], x[ok], (y + ci)[ok],
                        angle = 90, code = 3, length = 0.03, col = col, lwd = 1.6)
    lines(x, y, col = col, lwd = 2); points(x, y, col = col, pch = pch, bg = col, cex = 1.3)
  }
  eb(sub$recall_mean, sub$recall_ci, BLUE, 21)
  eb(sub$precision_mean, sub$precision_ci, ORANGE, 22)
  legend("bottomleft", c("recall", "precision"), col = c(BLUE, ORANGE),
         pch = c(21, 22), pt.bg = c(BLUE, ORANGE), lwd = 2, bty = "n", cex = 1.2)
}

roc_panel <- function(letter) {
  plot(NA, xlim = c(0, 1), ylim = c(0, 1.02), las = 1,
       xlab = "orphan false-positive rate", ylab = "planted true-host recall", main = "")
  title(main = letter, adj = 0, font.main = 2, cex.main = 1.7, line = 0.6)
  grid(col = "grey85", lty = 1)
  polygon(c(roc$fpr_grid, rev(roc$fpr_grid)), c(roc$ci_lo, rev(roc$ci_hi)),
          col = adjustcolor(BLUE, alpha.f = 0.18), border = NA)
  abline(0, 1, lty = 2, col = "grey60")
  lines(roc$fpr_grid, roc$mean_tpr, col = BLUE, lwd = 2.4)
  box()
  legend("bottomright", c("mean ROC", "95% CI across reps"),
         col = c(BLUE, adjustcolor(BLUE, alpha.f = 0.18)), lwd = c(2.4, 8),
         bty = "n", cex = 1.15)
  text(0.97, 0.34, sprintf("AUROC = %.2f", auroc), adj = 1, cex = 1.3, font = 2)
}

draw <- function() {
  par(mfrow = c(2, 2), mar = c(4.6, 4.6, 3, 1.4), oma = c(0, 0, 0.4, 0),
      cex.axis = 1.25, cex.lab = 1.4, mgp = c(2.9, 0.8, 0))
  metric_panel(d[d$axis == "complexity", ],    "community size (genomes)",   "a", logx = TRUE, shade = c(24, 58))
  metric_panel(d[d$axis == "coverage", ],       "donor sequencing depth (x)", "b")
  metric_panel(d[d$axis == "end2end_denovo", ], "community size (genomes)",   "c")
  roc_panel("d")
}

pdf(OUT, width = 11, height = 9.6); draw(); invisible(dev.off())
png(sub("\\.pdf$", ".png", OUT), width = 11, height = 9.6, units = "in", res = 150, type = "cairo")
draw(); invisible(dev.off())
cat("wrote", OUT, "| AUROC(mean curve) =", round(auroc, 3), "\n")
