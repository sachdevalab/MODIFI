#!/usr/bin/env Rscript
# Final combined figure (base R):
#   a  GC eta^2 vs flanking window (no zoom, no title)
#   b-d  context CV-R^2 heatmaps (C227, DSM, J99)
#   e-g  cross-species per-10-mer IPD scatter
#   h    metagenome cross-contig r by habitat
#   i    metagenome cross-contig r by taxonomic rank
# Run: conda run -n r_env Rscript benchmark/check_context/combined_final.R

CTX <- "/home/shuaiw/borg/revision/context"
FIG <- "/home/shuaiw/MODIFI/tmp/rev_figs"
out <- file.path(FIG, "combined_final.pdf")

col.C227 <- "#3b7dd8"; col.DSM <- "#2ca25f"; col.J99 <- "#d1495b"
pal.jet <- colorRampPalette(c("#00007F","#0000FF","#00A6FF","#27E36B",
                              "#B7F32C","#FFE31A","#F4E7C3"))(256)
sp.C227 <- bquote(italic("E. coli")~"C227")
sp.DSM  <- bquote(italic("C. israelensis")~"DSM 3043")
sp.J99  <- bquote(italic("H. pylori")~"J99")

addlab <- function(L) mtext(L, side=3, line=0.6, adj=-0.06, font=2, cex=1.4, xpd=NA)

## ---- readers -------------------------------------------------------------
gc <- read.csv(file.path(FIG, "gc_bin", "gc_bin_summary.csv"))

draw_gc <- function(){
  par(mar=c(4.2, 5.0, 1.6, 1.0))
  wins <- sort(unique(gc$flank_bp))
  plot(NA, xlim=range(log10(wins)), ylim=c(0, 100), axes=FALSE,
       xlab="", ylab="")
  for(s in list(c("E. coli C227",col.C227), c("C. israelensis DSM 3043",col.DSM),
                c("H. pylori J99",col.J99))){
    d <- gc[gc$strain==s[1],]; d <- d[order(d$flank_bp),]
    lines(log10(d$flank_bp), d$eta2_gc_pct, col=s[2], lwd=2)
    points(log10(d$flank_bp), d$eta2_gc_pct, col=s[2], pch=19, cex=1.1)
  }
  axis(1, at=log10(wins), labels=paste0("±", wins), cex.axis=0.95)
  axis(2, las=1, cex.axis=0.95); box()
  mtext("flanking window half-width (bp)", side=1, line=2.6, cex=0.95)
  mtext("within-10-mer IPD variation\nexplained by flanking GC  (%)", side=2, line=2.3, cex=0.9)
  legend("topleft", bty="n", cex=0.95, lwd=2, col=c(col.C227,col.DSM,col.J99),
         legend=as.expression(list(sp.C227, sp.DSM, sp.J99)))
  addlab("d")
}

draw_key <- function(){
  par(mar=c(3.0, 4.5, 2.6, 1.0))
  n <- length(pal.jet); xs <- seq(0.2, 0.8, length.out=n+1)
  plot.new(); plot.window(xlim=c(0.2, 0.2+3*(0.8-0.2)), ylim=c(0,1.3))
  rect(xs[-(n+1)], 0, xs[-1], 1, col=pal.jet, border=NA)
  rect(0.2,0,0.8,1, border="black")
  axis(1, at=c(0.2,0.4,0.6,0.8), line=-0.3, cex.axis=1.0)
  mtext(expression(R^2 ~ "value"), side=1, line=1.9, at=0.5, cex=0.95)
  text(0.2, 1.28, "Color Key", font=2, cex=1.25, adj=0, xpd=NA)
}

draw_heatmap <- function(strain_dir, titexpr, ylab, LET){
  m <- as.matrix(read.csv(file.path(CTX, strain_dir, "cv_r2_backoff_matrix.csv"),
                          row.names=1, check.names=FALSE))   # [d(6), u(11)]
  U <- 2:12; D <- 0:5
  mat <- t(m)                          # [u, d]
  z <- pmin(pmax(mat, 0.2), 0.8); z <- z[nrow(z):1,,drop=FALSE]
  vals <- mat[nrow(mat):1,,drop=FALSE]
  par(mar=c(3.2, if(ylab) 4.2 else 1.3, 2.4, 0.5))
  image(x=seq_along(D), y=seq_along(U), z=t(z), zlim=c(0.2,0.8), col=pal.jet,
        axes=FALSE, xlab="", ylab=""); box()
  for(i in seq_along(U)) for(j in seq_along(D)){
    v <- vals[i,j]
    text(j, i, sprintf("%.2f", v), cex=0.5, col=if(v<0.45) "white" else "black")
  }
  axis(1, at=seq_along(D), labels=paste0("-",D), tick=FALSE, line=-0.7, cex.axis=0.85)
  mtext("upstream context (5')", side=1, line=1.7, cex=0.75)
  if(ylab){
    axis(2, at=seq_along(U), labels=paste0("+",rev(U)), las=1, tick=FALSE,
         line=-0.4, cex.axis=0.85)
    mtext("downstream context (3')", side=2, line=2.7, cex=0.75)
  }
  title(main=titexpr, cex.main=1.15, font.main=1)
  addlab(LET)
}

draw_scatter <- function(f, xa, ya, xexpr, yexpr, LET){
  d <- read.csv(f)
  x <- d[[xa]]; y <- d[[ya]]
  r <- cor(x, y); n <- length(x)
  lim <- range(c(x,y)) + c(-0.2, 0.2)
  a <- min(0.55, max(0.03, 3500 / n))          # denser panels lighter
  cx <- if (n < 3000) 0.45 else 0.35           # sparse panels slightly bigger dots
  par(mar=c(4.2, 4.4, 2.0, 0.8))
  plot(x, y, pch=16, cex=cx, col=rgb(0,0,0,a), xlim=lim, ylim=lim,
       asp=1, xlab="", ylab="", las=1, cex.axis=0.9)
  abline(0, 1, col="#555", lwd=1)
  mtext(xexpr, side=1, line=2.5, cex=0.85)
  mtext(yexpr, side=2, line=2.4, cex=0.85)
  legend("topleft", bty="n", cex=0.95,
         legend=c(sprintf("r = %.2f", r), sprintf("n = %s", format(n, big.mark=","))))
  addlab(LET)
}

angle_axis <- function(labs, cex=0.85){
  u <- par("usr")
  text(seq_along(labs), u[3]-0.02*(u[4]-u[3]), labs, srt=30, adj=1, xpd=NA, cex=cex)
}

draw_box <- function(vals_list, labs, xlab, LET, ncnt){
  par(mar=c(5.5, 4.4, 2.0, 0.8))
  boxplot(vals_list, ylim=c(0,1), outline=FALSE, col="#9ecae1", border="#333",
          medcol="#d1495b", medlwd=2, whisklty=1, staplewex=0, xaxt="n",
          ylab="", las=1, cex.axis=0.95)
  mtext("Pearson r (per-10-mer IPD, contig pair)", side=2, line=2.6, cex=0.85)
  angle_axis(labs)
  u <- par("usr")
  text(seq_along(labs), 0.02, paste0("n=", ncnt), cex=0.7, col="#555")
  mtext(xlab, side=1, line=4.2, cex=0.9)
  addlab(LET)
}

## ---- assemble ------------------------------------------------------------
pdf(out, width=13, height=15.5)
par(oma=c(0, 0, 1.4, 0))          # room so top-row panel letters are not clipped
layout(matrix(c(1,1, 2,2,2,2,
                3,3, 4,4, 5,5,
                6,6, 7,7, 8,8,
                9,9,9, 10,10,10), nrow=4, byrow=TRUE),
       heights=c(0.8, 1.25, 1.15, 1.15))

draw_key()                                   # 1
draw_gc()                                    # 2
draw_heatmap("c224", sp.C227, TRUE,  "a")    # 3
draw_heatmap("dsm",  sp.DSM,  FALSE, "b")    # 4
draw_heatmap("j99",  sp.J99,  FALSE, "c")    # 5
draw_scatter(file.path(FIG,"cross_species_ipd","cross_species_C227_DSM.csv"),
             "z_C227","z_DSM", sp.C227, sp.DSM, "e")   # 6
draw_scatter(file.path(FIG,"cross_species_ipd","cross_species_C227_J99.csv"),
             "z_C227","z_J99", sp.C227, sp.J99, "f")   # 7
draw_scatter(file.path(FIG,"cross_species_ipd","cross_species_DSM_J99.csv"),
             "z_DSM","z_J99", sp.DSM, sp.J99, "g")     # 8

pairs <- read.csv(file.path(FIG,"meta_cross_contig","pairs.csv"))
# h: by habitat, ordered by median desc
med <- sort(tapply(pairs$pearson_r, pairs$environment, median), decreasing=TRUE)
hb <- names(med)
draw_box(lapply(hb, function(e) pairs$pearson_r[pairs$environment==e]), hb,
         "habitat", "h", sapply(hb, function(e) sum(pairs$environment==e)))
# i: by taxonomic rank
ro <- c("same_species","same_genus","same_family","same_order",
        "same_class","same_phylum","same_domain")
ro <- ro[ro %in% pairs$tax_rank]
dif <- c(same_genus="different_species", same_family="different_genus",
         same_order="different_family", same_class="different_order",
         same_phylum="different_class", same_domain="different_phylum")
tl <- sapply(ro, function(r) if(r %in% names(dif)) paste0(r,"\n(",dif[r],")") else r)
draw_box(lapply(ro, function(r) pairs$pearson_r[pairs$tax_rank==r]), tl,
         "lowest common taxonomic rank", "i",
         sapply(ro, function(r) sum(pairs$tax_rank==r)))

invisible(dev.off())
cat("wrote", out, "\n")
