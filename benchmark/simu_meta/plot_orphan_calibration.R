#!/usr/bin/env Rscript
# plot_orphan_calibration.R -- base-R threshold calibration figure (2x4), reading the Source
# Data CSVs so numbers match exactly. One panel per parameter: a final_score, b specificity,
# c min_sites, d motif fraction, e min_ctg_cov, f min_cov, g min_score. Each shows recall /
# precision / orphan false-positive rate (mean +/- 95% CI across the 5 orphan_300 replicates).
# Only panel letters are drawn; all other annotation lives in the caption. Points with no
# completed replicate (empty host_summary, e.g. min_score=100, min_ctg_cov=40) are dropped.

DIR  <- "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta/threshold"
OUT  <- file.path(DIR, "orphan_calibration.pdf")
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; RED <- "#CC3311"
rd <- function(f) read.csv(file.path(DIR, f), stringsAsFactors = FALSE)
sc  <- rd("orphan_calibration_score_sourcedata.csv")
sp  <- rd("orphan_calibration_specificity_sourcedata.csv")
ms  <- rd("orphan_calibration_minsites_sourcedata.csv")
mf  <- rd("orphan_calibration_minfrac_sourcedata.csv")
mcc <- rd("orphan_calibration_minctgcov_sourcedata.csv")
mcv <- rd("orphan_calibration_mincov_sourcedata.csv")
msc <- rd("orphan_calibration_minscore_sourcedata.csv")

letter <- function(l) title(main = l, adj = 0, font.main = 2, cex.main = 2.0, line = 0.6)
band <- function(x, m, ci, col) polygon(c(x, rev(x)), c(m - ci, rev(m + ci)),
                                        col = adjustcolor(col, 0.16), border = NA)

cont_panel <- function(d, xc, xlab, l, logx = FALSE) {
  x <- d[[xc]]
  plot(NA, xlim = range(x), ylim = c(0, 1.12), log = if (logx) "x" else "",
       xlab = xlab, ylab = "metric", las = 1)
  letter(l); grid(col = "grey88")
  band(x, d$recall_mean, d$recall_ci, BLUE)
  band(x, d$precision_mean, d$precision_ci, ORANGE)
  band(x, d$orphan_fpr_mean, d$orphan_fpr_ci, RED)
  lines(x, d$recall_mean, col = BLUE, lwd = 2.4)
  lines(x, d$precision_mean, col = ORANGE, lwd = 2.4)
  lines(x, d$orphan_fpr_mean, col = RED, lwd = 2.4, lty = 2)
  box()
}

disc_panel <- function(d, xc, xlab, l, logx = FALSE) {
  d <- d[!is.na(d$recall_mean), ]                 # drop empty (no completed rep) points
  x <- d[[xc]]
  plot(NA, xlim = range(x) + if (logx) 0 else c(-0.04, 0.04) * diff(range(x)),
       ylim = c(0, 1.12), log = if (logx) "x" else "",
       xlab = xlab, ylab = "metric", xaxt = "n", las = 1)
  letter(l); axis(1, at = x); grid(nx = NA, ny = NULL, col = "grey88")
  eb <- function(m, ci, col, pch, lty = 1) {
    ok <- ci > 0 & !is.na(ci)
    if (any(ok)) arrows(x[ok], (m - ci)[ok], x[ok], (m + ci)[ok], angle = 90, code = 3,
                        length = 0.03, col = col, lwd = 1.7)
    lines(x, m, col = col, lwd = 2.4, lty = lty); points(x, m, pch = pch, col = col, bg = col, cex = 1.5)
  }
  eb(d$recall_mean, d$recall_ci, BLUE, 21)
  eb(d$precision_mean, d$precision_ci, ORANGE, 22)
  eb(d$orphan_fpr_mean, d$orphan_fpr_ci, RED, 24, lty = 2)
  box()
}

draw <- function() {
  par(mfrow = c(2, 4), mar = c(4.8, 4.8, 3, 1.2), oma = c(0, 0, 0.4, 0),
      cex.axis = 1.2, cex.lab = 1.4, mgp = c(3.0, 0.9, 0))
  cont_panel(sc, "final_score", "final_score threshold (> x)", "a")
  legend(x = 0.02, y = 0.55, c("recall", "precision", "orphan-FPR"),
         col = c(BLUE, ORANGE, RED), lwd = 2.4, lty = c(1, 1, 2), bty = "n", cex = 1.15)
  cont_panel(sp, "specificity", "specificity threshold (< x)", "b", logx = TRUE)
  disc_panel(ms,  "min_sites",    "min modified-site count", "c")
  disc_panel(mf,  "min_frac",     "motif modified-fraction", "d")
  disc_panel(mcc, "min_ctg_cov",  "min contig coverage", "e")
  disc_panel(mcv, "min_cov",      "min per-site coverage", "f")
  disc_panel(msc, "min_score",    "min modification score", "g")
  plot.new()                                        # h: blank cell
}

pdf(OUT, width = 20, height = 10); draw(); invisible(dev.off())
png(sub("\\.pdf$", ".png", OUT), width = 20, height = 10, units = "in", res = 130, type = "cairo")
draw(); invisible(dev.off())
cat("wrote", OUT, "\n")
