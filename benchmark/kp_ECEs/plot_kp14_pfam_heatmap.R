#!/usr/bin/env Rscript
# Huge per-Pfam heatmap of the 14 Kp cluster representatives (rows) x every Pfam domain (columns).
# No side-annotation bars - clean matrix only. r_env Rscript.
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})
FIG <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
m <- as.matrix(read.csv(file.path(FIG,"kp14_pfam_matrix.csv"), row.names=1, check.names=FALSE))
title <- sprintf("Per-Pfam profile of the 14 Kp ECE representatives (%d Pfam domains)", ncol(m))
col_fun <- colorRamp2(c(0, 1, max(2,max(m))), c("#f7fbff","#9ecae1","#08306b"))
ht <- Heatmap(m, name="gene count", col=col_fun,
  cluster_rows=FALSE, cluster_columns=TRUE, show_column_dend=FALSE,
  row_names_side="left", row_names_gp=gpar(fontsize=9),
  column_names_gp=gpar(fontsize=2.4), column_names_rot=90,
  rect_gp=gpar(col="grey92", lwd=0.15),
  heatmap_legend_param=list(title="gene count", title_gp=gpar(fontsize=9,fontface="bold"),
                            labels_gp=gpar(fontsize=8)),
  width=unit(ncol(m)*0.115,"cm"), height=unit(nrow(m)*0.75,"cm"))
w <- ncol(m)*0.05 + 4        # inches; ~34 in wide (vector, zoomable)
h <- nrow(m)*0.32 + 3
drawf <- function(){ draw(ht, heatmap_legend_side="right", padding=unit(c(4,2,10,2),"mm"))
  grid.text(title, x=unit(0.5,"npc"), y=unit(1,"npc")-unit(4,"mm"), gp=gpar(fontsize=12,fontface="bold")) }
pdf(file.path(FIG,"fig_kp14_pfam_heatmap.pdf"), width=w, height=h); drawf(); dev.off()
png(file.path(FIG,"fig_kp14_pfam_heatmap.png"), width=w, height=h, units="in", res=150); drawf(); dev.off()
cat(sprintf("wrote fig_kp14_pfam_heatmap (%d x %d, %.0f x %.0f in)\n", nrow(m), ncol(m), w, h))
