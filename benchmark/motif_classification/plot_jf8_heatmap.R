#!/usr/bin/env Rscript
# JF8 (soil_1 control DB) modification-score heatmap: 14 real motifs x contigs,
# with a species color-annotation bar on the columns.
suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

datadir <- "/home/shuaiw/borg/revision/motif_class"
figdir  <- "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

mat <- as.matrix(read.csv(file.path(datadir, "jf8_heatmap_matrix.csv"),
                          row.names = 1, check.names = FALSE))
ann <- read.csv(file.path(datadir, "jf8_heatmap_colann.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
ann <- ann[match(colnames(mat), ann$contig), ]

# species colors (match the published panel)
sp_levels <- c("Ruminococcus gnavus", "Escherichia coli", "Clostridium bolteae",
               "Collinsella aerofaciens", "Bacteroides vulgatus",
               "Bacteroides thetaiotaomicron", "Bacteroides ovatus", "Bacteroides caccae")
sp_cols <- c("Ruminococcus gnavus"          = "#E41A1C",  # red
             "Escherichia coli"             = "#FFD700",  # gold
             "Clostridium bolteae"          = "#7FDB4A",  # light green
             "Collinsella aerofaciens"      = "#1B7837",  # dark green
             "Bacteroides vulgatus"         = "#33A1C9",  # teal
             "Bacteroides thetaiotaomicron" = "#5A9BD4",  # light blue
             "Bacteroides ovatus"           = "#7B3FA0",  # purple
             "Bacteroides caccae"           = "#000000")  # black

top_ann <- HeatmapAnnotation(
  Species = factor(ann$species, levels = sp_levels),
  col = list(Species = sp_cols),
  annotation_label = "Species",
  annotation_legend_param = list(
    Species = list(title = "Species", at = sp_levels, labels = sp_levels,
                   labels_gp = gpar(fontsize = 9, fontface = "italic"))),
  simple_anno_size = unit(4, "mm"),
  show_annotation_name = FALSE
)

col_fun <- colorRamp2(c(0, 0.5, 1), c("#f7f7f7", "#7f7f7f", "#000000"))

ht <- Heatmap(
  mat,
  name = "Modification\nscore (norm)",
  col = col_fun,
  top_annotation = top_ann,
  cluster_rows = TRUE, cluster_columns = TRUE,
  show_row_dend = TRUE, show_column_dend = TRUE,
  show_column_names = FALSE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),
  column_split = factor(ann$species, levels = sp_levels),
  column_title_gp = gpar(fontsize = 0),      # hide per-block titles (species shown by bar)
  column_gap = unit(0.6, "mm"),
  heatmap_legend_param = list(at = c(0, 0.5, 1)),
  border = TRUE
)

pdf(file.path(figdir, "jf8_heatmap_soil1.pdf"), width = 13.5, height = 5.2)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
dev.off()

png(file.path(figdir, "jf8_heatmap_soil1.png"), width = 13.5, height = 5.2, units = "in", res = 200)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
dev.off()

cat("wrote", file.path(figdir, "jf8_heatmap_soil1.pdf"), "and .png\n")
