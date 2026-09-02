#!/usr/bin/env Rscript
# plot_linkage_vs_completeness.R -- base-R figure (1x2) for ECE-host linkage across a
# gradient of assembled host quality, reading the Source Data TSVs written by
# plot_linkage_vs_completeness.py so numbers match exactly. Panel a: linkage rate vs. host
# completeness (all profiled host contigs, all 64 metagenomes). Panel b: linkage rate vs.
# host contamination (restricted to completeness >= 50% to avoid the completeness confound).
# Bars = fraction of host contigs linked to >=1 ECE; whiskers = Wilson 95% CI; n above each
# bar = contigs in that bin. Only panel letters are drawn; all other annotation lives in the
# caption.

DIR  <- "/home/shuaiw/MODIFI/tmp/rev_figs/complete_assess"
OUT  <- file.path(DIR, "linkage_vs_mag_quality.pdf")
BLUE <- "#0072B2"; GREY <- "#888888"

rd <- function(f) read.delim(file.path(DIR, f), stringsAsFactors = FALSE)
comp <- rd("linkage_vs_completeness_summary.tsv")
cont <- rd("linkage_vs_contamination_summary.tsv")

letter <- function(l) title(main = l, adj = 0, font.main = 2, cex.main = 2.0, line = 0.6)

rate_panel <- function(d, xlab, l) {
  ymax <- 0.25
  bp <- barplot(d$linkage_rate, names.arg = d$bin, col = BLUE, border = GREY,
                ylim = c(0, ymax), xlab = xlab,
                ylab = "Fraction of host contigs linked to >=1 ECE", las = 1)
  letter(l); box()
  arrows(bp, d$ci_low, bp, d$ci_high, angle = 90, code = 3, length = 0.04,
         col = "grey20", lwd = 1.7)
  text(bp, d$ci_high, labels = sprintf("%.3f\n(n=%d)", d$linkage_rate, d$n_contigs),
       pos = 3, offset = 0.4, cex = 1.05, xpd = NA)
}

draw <- function() {
  par(mfrow = c(1, 2), mar = c(4.8, 5, 3, 1.2), oma = c(0, 0, 0.4, 0),
      cex.axis = 1.25, cex.lab = 1.4, mgp = c(3.2, 0.9, 0))
  rate_panel(comp, "Host contig completeness (%)", "a")
  rate_panel(cont, "Host contig contamination (%)", "b")
}

pdf(OUT, width = 12, height = 5.2); draw(); invisible(dev.off())
png(sub("\\.pdf$", ".png", OUT), width = 12, height = 5.2, units = "in", res = 140,
    type = "cairo")
draw(); invisible(dev.off())
cat("wrote", OUT, "\n")
