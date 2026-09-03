#!/usr/bin/env Rscript
# Gene-profile heatmaps for the 14 Kp-cluster representatives: (1) curated ECE gene families,
# (2) COG functional categories (columns labelled with their meaning). Rows = representatives with
# side annotations (type, mobility, circular, known-ANI, AMR, toxin, virulence). r_env Rscript.
suppressMessages({library(ComplexHeatmap); library(circlize); library(grid)})
FIG <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
GP  <- "/home/shuaiw/borg/revision/kp_eces/gene_profile"

COG_NAMES <- c(J="Translation",A="RNA processing",K="Transcription",L="Replication/repair",
 B="Chromatin",D="Cell cycle/division",Y="Nuclear structure",V="Defense",T="Signal transduction",
 M="Cell wall/membrane",N="Cell motility",Z="Cytoskeleton",W="Extracellular",U="Secretion/trafficking",
 O="PTM/chaperones",X="Mobilome: phage/TE",C="Energy production",G="Carbohydrate",E="Amino acid",
 F="Nucleotide",H="Coenzyme",I="Lipid",P="Inorganic ion",Q="Secondary metabolites",
 R="General function",S="Function unknown")

ann <- read.delim(file.path(FIG,"kp14_annotation.tsv"), stringsAsFactors=FALSE)
ref <- read.delim(file.path(GP,"kp14_ani_known.tsv"), stringsAsFactors=FALSE)
ann <- merge(ann, ref[,c("representative","known_ANI")], by="representative", all.x=TRUE)
ann$rowlab <- sprintf("%s  (%.0f kb)", ann$cluster, ann$rep_length/1000)
rownames(ann) <- ann$representative
ann <- ann[order(-ann$rep_length),]

yn <- function(x) ifelse(x=="True"|x==TRUE|x=="yes","yes","no")
legpar <- list(title_gp=gpar(fontsize=9,fontface="bold"), labels_gp=gpar(fontsize=8),
               grid_height=unit(3.5,"mm"), grid_width=unit(3.5,"mm"))

draw_hm <- function(mat_csv, outbase, title, cluster_cols, col_label_map=NULL){
  m <- as.matrix(read.csv(mat_csv, row.names=1, check.names=FALSE))[ann$representative,,drop=FALSE]
  rownames(m) <- ann$rowlab
  collab <- colnames(m)
  if(!is.null(col_label_map)) collab <- paste0(colnames(m), ": ", col_label_map[colnames(m)])
  mob <- ann$predicted_mobility; mob[is.na(mob)|mob==""] <- "unknown"
  namr <- suppressWarnings(as.integer(ann$n_amr)); namr[is.na(namr)] <- 0
  la <- rowAnnotation(
    Type = ann$rep_type,
    Mobility = mob,
    Circular = yn(ann$rep_circular),
    `Known (ANI)` = yn(ann$known_ANI),
    AMR = ifelse(namr>0,"yes","no"),
    Toxin = yn(ann$has_toxin),
    Virulence = yn(ann$has_virulence),
    col = list(Type=c(plasmid="#2171b5", virus="#31a354"),
               Mobility=c(conjugative="#cb181d", mobilizable="#fd8d3c", `non-mobilizable`="#bdbdbd", unknown="#f0f0f0"),
               Circular=c(yes="#525252", no="#e0e0e0"),
               `Known (ANI)`=c(yes="#6a51a3", no="#e0e0e0"),
               AMR=c(yes="#b8860b", no="#f0f0f0"),
               Toxin=c(yes="#238b8b", no="#f0f0f0"),
               Virulence=c(yes="#ce1256", no="#f0f0f0")),
    annotation_name_gp=gpar(fontsize=8.5), annotation_name_side="top",
    simple_anno_size=unit(4,"mm"), gap=unit(0.6,"mm"),
    annotation_legend_param=list(Type=legpar,Mobility=legpar,Circular=legpar,
                                 `Known (ANI)`=legpar,AMR=legpar,Toxin=legpar,Virulence=legpar))
  ra <- rowAnnotation(`total genes`=anno_barplot(rowSums(m), gp=gpar(fill="#6baed6"),
                        width=unit(2,"cm"), axis_param=list(gp=gpar(fontsize=7))),
                        annotation_name_side="bottom")
  col_fun <- colorRamp2(c(0, max(1,max(m))/2, max(1,max(m))), c("#f7fbff","#6baed6","#08306b"))
  ht <- Heatmap(m, name="gene count", col=col_fun, column_labels=collab,
    left_annotation=la, right_annotation=ra,
    cluster_rows=FALSE, cluster_columns=cluster_cols,
    row_names_side="left", row_names_gp=gpar(fontsize=9),
    column_names_gp=gpar(fontsize=9), column_names_rot=45,
    cell_fun=function(j,i,x,y,w,h,fill){v<-m[i,j]; if(v>0) grid.text(v,x,y,gp=gpar(fontsize=7,col=ifelse(v>max(m)*0.6,"white","grey20")))},
    column_title=NULL,
    heatmap_legend_param=c(list(title="gene count"),legpar),
    width=unit(ncol(m)*0.75,"cm"), height=unit(nrow(m)*0.75,"cm"))
  extra <- if(!is.null(col_label_map)) 4 else 0
  w <- 8 + ncol(m)*0.5 + extra; h <- 5 + nrow(m)*0.42
  pdfpath <- file.path(FIG,paste0(outbase,".pdf"))
  pdf(pdfpath, width=w, height=h)
  draw(ht, merge_legend=TRUE, heatmap_legend_side="right", annotation_legend_side="right",
       padding=unit(c(4,2,16,2),"mm")); dev.off()
  # crop the oversized canvas to the content bbox (+ matching PNG)
  system(paste("bash /home/shuaiw/MODIFI/benchmark/kp_ECEs/crop_pdf.sh", shQuote(pdfpath)))
  write.csv(m, file.path(FIG,paste0(outbase,"_sourcedata.csv")))
  cat("wrote",outbase,"\n")
}

draw_hm(file.path(FIG,"kp14_family_matrix.csv"), "fig_kp14_gene_profile_families",
        "Gene-family profile: 14 Kp ECE representatives", FALSE)
draw_hm(file.path(FIG,"kp14_cog_matrix.csv"), "fig_kp14_gene_profile_cog",
        "COG-category profile: 14 Kp ECE representatives", TRUE, COG_NAMES)
