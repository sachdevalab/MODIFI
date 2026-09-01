#!/usr/bin/env Rscript
# plot_strain_benchmark.R -- combined strain-level benchmark figure (base R), 2 rows x 3 cols.
# Merges the strain-mixture (C2) panels and the E. coli strain-mixture assembly panel. Reads the
# Source Data CSVs so numbers match the Python figures. Only panel letters are drawn; titles and
# value labels are removed (they live in the caption). Species names are italic.

BASE <- "/home/shuaiw/MODIFI/tmp/rev_figs/simu_meta"
SM   <- read.csv(file.path(BASE, "strain_het", "strain_mix_recall_precision_sourcedata.csv"),
                 stringsAsFactors = FALSE)
EC   <- read.csv(file.path(BASE, "e_coli", "ecoli_panel_sourcedata.csv"), stringsAsFactors = FALSE)
MR   <- read.csv(file.path(BASE, "strain_het", "strain_mix_motif_recovery_sourcedata.csv"),
                 stringsAsFactors = FALSE)
OUT  <- file.path(BASE, "strain_het", "strain_benchmark_combined.pdf")
BLUE <- "#0072B2"; ORANGE <- "#D55E00"; GREEN <- "#009E73"; GREY <- "#BBBBBB"

ecv <- setNames(EC$value, EC$metric)
klab <- paste0("K=", SM$K)
x <- seq_len(nrow(SM))
XL <- "strains per donor species (K)"

letter <- function(l) title(main = l, adj = 0, font.main = 2, cex.main = 2.1, line = 0.7)

line_panel <- function(rm, rc, pm, pc, ylab, l, ylim = c(0, 1.06)) {
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = ylim, xaxt = "n",
       xlab = XL, ylab = ylab, las = 1)
  letter(l); axis(1, at = x, labels = klab)
  grid(nx = NA, ny = NULL, col = "grey85")
  eb <- function(y, ci, col, pch) {
    ok <- ci > 0
    if (any(ok)) arrows(x[ok], (y - ci)[ok], x[ok], (y + ci)[ok],
                        angle = 90, code = 3, length = 0.03, col = col, lwd = 1.7)
    lines(x, y, col = col, lwd = 2.2); points(x, y, pch = pch, col = col, bg = col, cex = 1.5)
  }
  eb(rm, rc, BLUE, 21); eb(pm, pc, ORANGE, 22)
  legend("bottomleft", c("recall", "precision"), col = c(BLUE, ORANGE),
         pch = c(21, 22), pt.bg = c(BLUE, ORANGE), lwd = 2, bty = "n", cex = 1.35)
}

draw <- function() {
  par(mfrow = c(2, 4), mar = c(6, 5, 3, 1.2), oma = c(0, 0, 0.4, 0),
      cex.axis = 1.3, cex.lab = 1.5, mgp = c(3.1, 0.9, 0))

  # a: community strain composition (grouped bars)
  comp <- rbind(SM$donor_strains, SM$total_genomes)
  barplot(comp, beside = TRUE, col = c(BLUE, GREY), names.arg = klab, border = NA,
          ylab = "count", ylim = c(0, max(comp) * 1.12), las = 1)
  letter("a"); box()
  legend("topleft", c("donor strains", "total genomes (incl. 242 bg)"),
         fill = c(BLUE, GREY), border = NA, bty = "n", cex = 1.3)

  # b: motif-detection recovery vs K (isolate motifs = ground truth; same filter both sides)
  xm <- seq_len(nrow(MR))
  plot(NA, xlim = c(0.8, length(xm) + 0.2), ylim = c(0, 1.06), xaxt = "n",
       xlab = XL, ylab = "motif recovery rate", las = 1)
  letter("b"); axis(1, at = xm, labels = paste0("K=", MR$K)); grid(nx = NA, ny = NULL, col = "grey85")
  okm <- MR$recovery_pooled_ci > 0
  arrows(xm[okm], (MR$recovery_pooled_mean - MR$recovery_pooled_ci)[okm], xm[okm],
         (MR$recovery_pooled_mean + MR$recovery_pooled_ci)[okm], angle = 90, code = 3,
         length = 0.03, col = GREEN, lwd = 1.7)
  lines(xm, MR$recovery_pooled_mean, col = GREEN, lwd = 2.2)
  points(xm, MR$recovery_pooled_mean, pch = 21, col = GREEN, bg = GREEN, cex = 1.5)

  # c: species-level recall & precision
  line_panel(SM$species_recall_mean, SM$species_recall_ci,
             SM$species_precision_mean, SM$species_precision_ci, "species-level metric", "c")

  # d: strain-level recall & precision
  line_panel(SM$strain_recall_mean, SM$strain_recall_ci,
             SM$strain_precision_mean, SM$strain_precision_ci, "strain-level metric", "d")

  # e: strain accuracy (%)
  plot(NA, xlim = c(0.8, length(x) + 0.2), ylim = c(80, 100.5), xaxt = "n",
       xlab = XL, ylab = "strain accuracy (%)", las = 1)
  letter("e"); axis(1, at = x, labels = klab); grid(nx = NA, ny = NULL, col = "grey85")
  ok <- SM$strain_accuracy_ci > 0
  arrows(x[ok], (SM$strain_accuracy_mean - SM$strain_accuracy_ci)[ok], x[ok],
         (SM$strain_accuracy_mean + SM$strain_accuracy_ci)[ok], angle = 90, code = 3,
         length = 0.03, col = GREEN, lwd = 1.7)
  lines(x, SM$strain_accuracy_mean, col = GREEN, lwd = 2.2)
  points(x, SM$strain_accuracy_mean, pch = 21, col = GREEN, bg = GREEN, cex = 1.5)

  # f: E. coli strain mixture -> chimeric assembly (italic species name, 2-line labels below axis)
  bp <- barplot(c(ecv["chimeric_ecoli_contigs_pct"], ecv["checkm2_completeness_pct"]),
                col = c(ORANGE, GREY), names.arg = c("", ""), border = NA,
                ylab = "percent", ylim = c(0, 108), las = 1, width = 0.7, space = 0.8)
  letter("f"); box()
  yb <- -108 * 0.055
  text(bp[1], yb, expression(atop("chimeric", italic("E. coli") ~ "contigs")), xpd = NA, cex = 1.3)
  text(bp[2], yb, expression(atop("CheckM2", "completeness")), xpd = NA, cex = 1.3)

  # g: E. coli de-novo ECE->strain linkage
  barplot(c(ecv["denovo_strain_recall"], ecv["denovo_strain_precision"]),
          col = c(BLUE, ORANGE), names.arg = c("recall", "precision"), border = NA,
          ylab = "strain-level metric", ylim = c(0, 1.1), las = 1, width = 0.7, space = 0.8)
  letter("g"); box()
  plot.new()                                              # h: blank cell
}

pdf(OUT, width = 22, height = 10); draw(); invisible(dev.off())
png(sub("\\.pdf$", ".png", OUT), width = 22, height = 10, units = "in", res = 130, type = "cairo")
draw(); invisible(dev.off())
cat("wrote", OUT, "\n")
