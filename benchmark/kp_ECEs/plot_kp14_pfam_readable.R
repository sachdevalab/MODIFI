#!/usr/bin/env Rscript
# ONE combined shared-Pfam figure, two panels in a single ROW (side by side): (left) curated
# functional families with a left family colour strip, split by family; (right) Other/hypothetical
# as its own panel. Both panels share the same 14-rep column order + plasmid/virus Type bar; a single
# shared legend block (ECE type / functional family / gene count) sits at the far right.
# Rows = Pfams (>=5 reps), cols = 14 reps. No dendrogram.
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})
FIG <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
GP  <- "/home/shuaiw/borg/revision/kp_eces/gene_profile"
m <- as.matrix(read.csv(file.path(FIG,"kp14_pfam_matrix.csv"), row.names=1, check.names=FALSE))
fam <- read.csv(file.path(FIG,"kp14_pfam_family_map.csv"), stringsAsFactors=FALSE)
famv <- setNames(fam$family, fam$pfam)
m <- m[, colnames(m)[colSums(m>0) >= 5], drop=FALSE]
mt <- t(m)                                              # Pfam x rep

rep <- read.csv(file.path(GP,"kp14_representatives.tsv"), sep="\t", stringsAsFactors=FALSE)
rep$lab <- sprintf("%s  (%.0f kb)", rep$cluster, rep$rep_length/1000)
coltype <- setNames(rep$rep_type, rep$lab)[colnames(mt)]; coltype[is.na(coltype)] <- "plasmid"
colord <- colnames(mt)[hclust(dist(t(mt)))$order]      # shared column order for both panels

FAM_COL <- c("Replication"="#1f78b4","Partitioning/stability"="#a6cee3","Conjugation/mobilization"="#e31a1c",
 "Transposon/integrase"="#6a3d9a","Defense (RM/CRISPR/Abi)"="#ff7f00","Toxin-antitoxin"="#238b8b",
 "Metal/biocide resistance"="#b15928","AMR (antibiotic)"="#fb9a99","Virulence"="#ce1256",
 "Phage structural"="#33a02c","Metabolism"="#b2df8a","Other/hypothetical"="#d9d9d9")
legpar <- list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8))
funcorder <- c("Replication","Partitioning/stability","Conjugation/mobilization","Transposon/integrase",
 "Defense (RM/CRISPR/Abi)","Toxin-antitoxin","Metal/biocide resistance","AMR (antibiotic)","Virulence",
 "Phage structural","Metabolism")

func_rows  <- rownames(mt)[famv[rownames(mt)] != "Other/hypothetical"]
other_rows <- rownames(mt)[famv[rownames(mt)] == "Other/hypothetical"]
mx <- max(mt[c(func_rows,other_rows),])                 # shared colour scale
col_fun <- colorRamp2(c(0,1,max(2,mx)), c("#f7fbff","#9ecae1","#08306b"))
mkcell <- function(x) function(j,i,xx,yy,w,hh,f){v<-x[i,j]; if(v>0) grid.text(v,xx,yy,gp=gpar(fontsize=4.2,col=ifelse(v>mx*0.6,"white","grey25")))}
mktop  <- function(x) HeatmapAnnotation(Type=coltype[colnames(x)],
             col=list(Type=c(plasmid="#2171b5", virus="#31a354")),
             annotation_name_gp=gpar(fontsize=8), annotation_name_side="right", show_legend=FALSE)

build_ht <- function(rows, split_by_family){
  x <- mt[rows,,drop=FALSE]
  args <- list(x, name=if(split_by_family)"gene count" else "gc2", col=col_fun,
    column_order=colord, cluster_columns=FALSE, cluster_rows=TRUE, show_row_dend=FALSE,
    show_heatmap_legend=FALSE, top_annotation=mktop(x),
    row_names_gp=gpar(fontsize=5.6),
    column_names_side="bottom", column_names_gp=gpar(fontsize=7), column_names_rot=45,
    rect_gp=gpar(col="grey90", lwd=0.2), cell_fun=mkcell(x))
  if(split_by_family){
    args <- c(args, list(row_names_side="left",
      left_annotation=rowAnnotation(Family=famv[rows], col=list(Family=FAM_COL),
        show_annotation_name=FALSE, show_legend=FALSE, width=unit(3,"mm")),
      row_split=factor(famv[rows], levels=intersect(funcorder, unique(famv[rows]))),
      cluster_row_slices=FALSE, row_gap=unit(0.8,"mm"), row_title=NULL))
  } else {
    args <- c(args, list(row_names_side="right",
      row_split=factor(rep("Other/hypothetical", length(rows))),
      row_title="Other / hypothetical", row_title_rot=90,
      row_title_gp=gpar(fontsize=8,fontface="bold"), row_title_side="left"))
  }
  do.call(Heatmap, args)
}

ht_func  <- build_ht(func_rows,  TRUE)
ht_other <- build_ht(other_rows, FALSE)

fam_used <- intersect(funcorder, unique(famv[func_rows]))
lg_type <- Legend(title="ECE type", labels=c("plasmid","virus"),
                  legend_gp=gpar(fill=c("#2171b5","#31a354")),
                  title_gp=legpar$title_gp, labels_gp=legpar$labels_gp)
lg_fam  <- Legend(title="functional family", labels=fam_used, legend_gp=gpar(fill=FAM_COL[fam_used]),
                  title_gp=legpar$title_gp, labels_gp=legpar$labels_gp)
lg_gc   <- Legend(title="gene count", col_fun=col_fun, at=pretty(c(0,mx)),
                  title_gp=legpar$title_gp, labels_gp=legpar$labels_gp)
lgs <- packLegend(lg_type, lg_fam, lg_gc, gap=unit(6,"mm"))

outbase <- "fig_kp14_pfam_shared"; W <- 13; H <- 8.6
for(dev in c("pdf","png")){
  if(dev=="pdf") pdf(file.path(FIG,paste0(outbase,".pdf")), width=W, height=H)
  else png(file.path(FIG,paste0(outbase,".png")), width=W, height=H, units="in", res=200)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,3,
      widths=unit.c(unit(0.44,"npc"),unit(0.40,"npc"),unit(0.16,"npc")))))
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=1)); draw(ht_func,  newpage=FALSE); popViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=2)); draw(ht_other, newpage=FALSE); popViewport()
  pushViewport(viewport(layout.pos.row=1, layout.pos.col=3)); draw(lgs, x=unit(2,"mm"), just="left"); popViewport()
  popViewport(); dev.off()
}
write.csv(mt[c(func_rows,other_rows),], file.path(FIG,paste0(outbase,"_sourcedata.csv")))
cat("wrote",outbase,": functional",length(func_rows),"+ other",length(other_rows),"Pfams\n")
