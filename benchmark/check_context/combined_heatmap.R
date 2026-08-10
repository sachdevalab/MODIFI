#!/usr/bin/env Rscript
# Combined 3-strain k-mer context R^2 figure with one shared color key.
# Reads the back-off CV-R^2 matrices (rows d=0..5, cols u=2..12) for C227, DSM 3043
# and J99, and draws three heatmaps in the reference style:
#   y-axis = downstream context (3'), labelled +2..+12 bp (top -> bottom)
#   x-axis = upstream   context (5'), labelled -0..-5 bp (left -> right)
# Shared colour scale 0.2-0.8 (jet -> pale tan). PDF only.
#
# Run:  conda run -n r_env Rscript benchmark/check_context/combined_heatmap.R

## ---- inputs ---------------------------------------------------------------
base <- "/home/shuaiw/borg/revision/context"
strains <- list(
  list(file = file.path(base, "c224", "cv_r2_backoff_matrix.csv"),
       title = "E. coli C227"),
  list(file = file.path(base, "dsm", "cv_r2_backoff_matrix.csv"),
       title = "C. israelensis DSM 3043"),
  list(file = file.path(base, "j99", "cv_r2_backoff_matrix.csv"),
       title = "H. pylori J99")
)
out_pdf <- "/home/shuaiw/MODIFI/tmp/rev_figs/combined_context_r2.pdf"

ZMIN <- 0.2; ZMAX <- 0.8
U <- 2:12          # upstream extent  -> rows, "+N bp"
D <- 0:5           # downstream extent -> cols, "-N bp"

## ---- palette: deep blue -> blue -> cyan -> green -> yellow -> pale tan -----
pal <- colorRampPalette(c("#00007F", "#0000FF", "#00A6FF", "#27E36B",
                          "#B7F32C", "#FFE31A", "#F4E7C3"))(256)

read_mat <- function(f) {
  # CSV is [d rows, u cols]; return display matrix [u, d] (rows=u, cols=d)
  m <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE))  # 6 x 11 (d x u)
  t(m)                                                             # 11 x 6 (u x d)
}

draw_panel <- function(mat, title, show_ylab) {
  # mat: [u (11), d (6)]. image() needs z[x, y]; we want x=d (cols), y=u (rows),
  # with +2 (u=2) at TOP. Put u on y and reverse so small u is at the top.
  z <- pmin(pmax(mat, ZMIN), ZMAX)          # clamp for colouring
  z <- z[nrow(z):1, , drop = FALSE]         # reverse u so +2 ends up on top
  vz <- mat[nrow(mat):1, , drop = FALSE]    # unclamped values, same orientation
  zt <- t(z)                                # image z[x=d, y=u]
  par(mar = c(3.6, if (show_ylab) 4.6 else 1.0, 2.6, 0.6))
  image(x = seq_along(D), y = seq_along(U), z = zt,
        zlim = c(ZMIN, ZMAX), col = pal, axes = FALSE, xlab = "", ylab = "")
  # per-cell value labels (white ink on the dark low-R^2 cells, black elsewhere)
  for (p in seq_along(U)) for (j in seq_along(D))
    text(j, p, sprintf("%.3f", vz[p, j]), cex = 0.78,
         col = if (vz[p, j] < 0.48) "white" else "black")
  box(lwd = 1)
  # x ticks: -0 .. -5 (upstream 5')
  axis(1, at = seq_along(D), labels = paste0("-", D), tick = FALSE,
       line = -0.4, cex.axis = 1.05)
  mtext("upstream context (5')", side = 1, line = 2.3, cex = 0.95)
  # y ticks: +2 (top) .. +12 (bottom)  -> after the reverse, top row = u=2
  if (show_ylab) {
    axis(2, at = seq_along(U), labels = paste0("+", rev(U)), las = 1,
         tick = FALSE, line = -0.3, cex.axis = 1.05)
    mtext("downstream context (3')", side = 2, line = 3.1, cex = 0.95)
  }
  title(main = title, font.main = 3, cex.main = 1.35)
}

draw_key <- function() {
  # horizontal colour key: "Color Key" title, bar, ticks, R^2 label; left ~1/3 width
  par(mar = c(3.2, 5.2, 2.8, 1.0))
  n <- length(pal)
  xs <- seq(ZMIN, ZMAX, length.out = n + 1)
  plot.new(); plot.window(xlim = c(ZMIN, ZMIN + 3 * (ZMAX - ZMIN)), ylim = c(0, 1.3))
  rect(xs[-(n + 1)], 0, xs[-1], 1, col = pal, border = NA)
  rect(ZMIN, 0, ZMAX, 1, border = "black", lwd = 1)
  ticks <- c(0.2, 0.4, 0.6, 0.8)
  axis(1, at = ticks, labels = ticks, line = -0.2, cex.axis = 1.15)
  mtext(expression(R^2 ~ "value"), side = 1, line = 2.1, at = mean(c(ZMIN, ZMAX)),
        cex = 1.1)
  text(ZMIN, 1.28, "Color Key", xpd = NA, font = 2, cex = 1.5, adj = 0)
}

## ---- assemble -------------------------------------------------------------
mats <- lapply(strains, function(s) read_mat(s$file))

pdf(out_pdf, width = 13, height = 6.6)
layout(matrix(c(1, 1, 1,
                2, 3, 4), nrow = 2, byrow = TRUE), heights = c(0.26, 1))
draw_key()
for (i in seq_along(strains))
  draw_panel(mats[[i]], strains[[i]]$title, show_ylab = (i == 1))
invisible(dev.off())

cat("wrote", out_pdf, "\n")
for (i in seq_along(strains))
  cat(sprintf("  %-26s plateau(max)=%.3f  [-2,+7]=%.3f\n",
              strains[[i]]$title, max(mats[[i]]),
              mats[[i]][which(U == 7), which(D == 2)]))
