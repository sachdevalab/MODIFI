#!/usr/bin/env Rscript

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)

plot_heatmap <- function(data_file, heat_map) {
  # Read updated data file (has sample, type, anno, contig, and updated motif fractions)
  df <- read.csv(data_file)
  
  cat("Data shape:", nrow(df), "x", ncol(df), "\n")
  cat("Columns:", paste(colnames(df), collapse=", "), "\n")
  print(head(df))
  
  # Get motif columns (exclude metadata columns)
  motif_columns <- setdiff(colnames(df), c('sample', 'type', 'anno', 'contig'))
  cat("Motif columns:", length(motif_columns), "\n")
  
  # Sort data by sample and type order (host=0, P1=1, P2=2)
  df_sorted <- df %>%
    mutate(
      # Fix anno: change P2 to P1 where type is plasmid_1
      anno = ifelse(grepl('plasmid_1', type) & grepl('_P2$', anno), 
                    sub('_P2$', '_P1', anno), 
                    anno),
      type_order = case_when(
        grepl('host', type) ~ 0,
        grepl('P1', anno) ~ 1,
        grepl('P2', anno) ~ 2,
        TRUE ~ 3
      )
    ) %>%
    arrange(sample, type_order) %>%
    select(-type_order)
  
  # Create unique row names using anno (make unique if duplicates exist)
  row_labels <- make.unique(as.character(df_sorted$anno), sep = "_")
  rownames(df_sorted) <- row_labels
  
  # Extract heatmap data matrix
  heatmap_data <- as.matrix(df_sorted[, motif_columns])
  rownames(heatmap_data) <- row_labels
  
  # Create type and sample labels
  type_labels <- sapply(df_sorted$type, function(type_str) {
    if (grepl('plasmid', type_str)) return('plasmid')
    else if (grepl('host', type_str)) return('host')
    else return('unknown')
  })
  sample_labels <- df_sorted$sample
  
  type_colors <- c('plasmid' = '#66c2a5', 'host' = '#fc8d62', 'unknown' = '#bbbbbb')
  
  # Create better sample colors using viridis-like palette
  unique_samples <- unique(sample_labels)
  n_samples <- length(unique_samples)
  
  if (n_samples <= 6) {
    # Use distinct colors for few samples
    sample_color_vals <- c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02')[1:n_samples]
  } else {
    # Use Set3 for more samples
    sample_color_vals <- brewer.pal(min(max(3, n_samples), 12), "Set3")[1:n_samples]
  }
  
  sample_colors <- setNames(sample_color_vals, unique_samples)
  
  # Create row annotation (only Type; sample annotation removed)
  row_ha <- rowAnnotation(
    Type = type_labels,
    col = list(
      Type = type_colors
    ),
    show_annotation_name = FALSE,
    annotation_name_side = "top",
    annotation_legend_param = list(
      Type = list(
        title = "",
        ncol = 2,
        title_gp = gpar(fontsize = 0),
        labels_gp = gpar(fontsize = 7),
        grid_height = unit(3, "mm"),
        grid_width = unit(3, "mm")
      )
    )
  )
  
  # Create heatmap with square cells
  cell_size <- unit(5, "mm")
  heatmap_args <- list(
    matrix = heatmap_data,
    name = "Modification Fraction",
    col = colorRamp2(
      c(0, 0.5, 1),
      c("white", "#fc8d59", "#67001f")
    ),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    clustering_method_columns = "average",
    clustering_distance_columns = "euclidean",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    column_names_rot = 45,
    column_names_side = "bottom",
    left_annotation = row_ha,
    width = ncol(heatmap_data) * cell_size,
    height = nrow(heatmap_data) * cell_size,
    column_dend_height = unit(4, "mm"),
    column_names_max_height = unit(3, "cm"),
    row_names_max_width = unit(5, "cm"),
    heatmap_legend_param = list(
      title = "Modification Fraction",
      direction = "horizontal",
      legend_width = unit(2, "cm"),
      grid_height = unit(3, "mm"),
      title_gp = gpar(fontsize = 6),
      labels_gp = gpar(fontsize = 6)
    ),
    rect_gp = gpar(col = "gray90", lwd = 0.5)
  )
  
  ht <- do.call(Heatmap, heatmap_args)
  
  # Save to file
  pdf(heat_map, width = 6, height = 5)
  draw(ht, heatmap_legend_side = "top", annotation_legend_side = "top")
  dev.off()
  
  cat(sprintf("\nHeatmap saved to: %s\n", heat_map))
}

# Main execution
data_file <- "180_4_updated.csv"
heat_map <- "../../tmp/figures/inversion/motif_change_heatmap.pdf"

plot_heatmap(data_file, heat_map)
