#!/usr/bin/env Rscript

library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(tidyr)

heatmap_plot <- function(df, heat_map, motifs_to_remove = c()) {
  # Check if data is in long format (contig, motifString, fraction)
  if ("motifString" %in% colnames(df) && "fraction" %in% colnames(df)) {
    # Filter out unwanted motifs
    if (length(motifs_to_remove) > 0) {
      df <- df %>%
        filter(!motifString %in% motifs_to_remove)
      cat(sprintf("Removed %d motifs: %s\n", 
                  length(motifs_to_remove), 
                  paste(motifs_to_remove, collapse=", ")))
    }
    
    # Pivot to wide format: rows = contigs, columns = motifs
    df_wide <- df %>%
      select(contig, motifString, fraction) %>%
      pivot_wider(names_from = motifString, values_from = fraction, values_fill = 0)
    
    heatmap_data <- df_wide %>%
      select(-contig) %>%
      as.matrix()
    
    rownames(heatmap_data) <- df_wide$contig
    
    # Create heatmap without annotations (simple version)
    ht <- Heatmap(
      heatmap_data,
      name = "Fraction",
      col = colorRamp2(c(0, 0.3, 0.7, 1), c("white", "#fee5d9", "#de2d26", "#67001f")),
      
      # Clustering
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      clustering_distance_rows = "euclidean",
      clustering_distance_columns = "euclidean",
      clustering_method_rows = "average",
      clustering_method_columns = "average",
      
      # Cell borders
      border = TRUE,
      rect_gp = gpar(col = "gray90", lwd = 0.5),
      
      # Row/column names
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_side = "right",
      row_names_gp = gpar(fontsize = 9),
      column_names_gp = gpar(fontsize = 7),
      column_names_rot = 90,
      column_names_centered = FALSE,
      
      # Size
      width = unit(14, "cm"),
      height = unit(0.4 * nrow(heatmap_data), "cm"),
      
      # Legend
      heatmap_legend_param = list(
        title = "Fraction",
        title_gp = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        direction = "vertical",
        legend_height = unit(4, "cm")
      )
    )
    
  } else {
    # Original code for wide format with sample/type/anno columns
    motif_columns <- setdiff(colnames(df), c('sample', 'type', 'anno', 'contig'))
    
    # Sort data by sample and type order
    df_sorted <- df %>%
      mutate(
        type_order = case_when(
          grepl('host', type) ~ 0,
          grepl('P1', anno) ~ 1,
          grepl('P2', anno) ~ 2,
          TRUE ~ 3
        )
      ) %>%
      arrange(sample, type_order) %>%
      select(-type_order)
    
    # Prepare heatmap data matrix
    heatmap_data <- df_sorted %>%
      select(all_of(motif_columns)) %>%
      as.matrix()
    
    rownames(heatmap_data) <- df_sorted$anno
  
  # Create type colors
  get_base_type <- function(type_str) {
    if (grepl('plasmid', type_str)) {
      return('plasmid')
    } else if (grepl('host', type_str)) {
      return('host')
    } else {
      return('unknown')
    }
  }
  
  type_labels <- sapply(df_sorted$type, get_base_type)
  sample_labels <- df_sorted$sample
  
  # Define color palettes
  type_palette <- c('plasmid' = '#1f77b4', 'host' = '#ff7f0e', 'unknown' = '#bbbbbb')
  type_colors <- type_palette[type_labels]
  
  # Sample colors using Set3-like palette
  unique_samples <- unique(sample_labels)
  sample_palette <- hcl.colors(length(unique_samples), palette = "Set3")
  names(sample_palette) <- unique_samples
  sample_colors <- sample_palette[sample_labels]
  
  # Create row annotations
  row_anno <- rowAnnotation(
    Type = type_labels,
    Sample = sample_labels,
    col = list(
      Type = type_palette,
      Sample = sample_palette
    ),
    annotation_width = unit(c(3, 3), "mm"),
    show_annotation_name = TRUE,
    show_legend = c(TRUE, FALSE),
    annotation_name_side = "top"
  )
  
  # Create heatmap
  ht <- Heatmap(
    heatmap_data,
    name = "Value",
    col = colorRamp2(seq(min(heatmap_data), max(heatmap_data), length = 100), 
                     c("#440154FF", "#31688EFF", "#35B779FF", "#FDE724FF")),
    
    # Clustering
    cluster_rows = FALSE,  # Don't cluster rows to keep sample grouping
    cluster_columns = TRUE,  # Cluster columns (motifs)
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "average",
    
    # Row/column names
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    
    # Row annotation
    left_annotation = row_anno,
    
    # Size
    width = unit(12, "cm"),
    height = unit(10, "cm"),
    
    # Legend
    heatmap_legend_param = list(
      title = "Value",
      direction = "vertical"
    )
  )
  }
  
  # Save to file - dynamic sizing based on data dimensions
  if ("motifString" %in% colnames(df) && "fraction" %in% colnames(df)) {
    pdf_width <- min(max(8, ncol(heatmap_data) * 0.15), 20)
    pdf_height <- min(max(6, nrow(heatmap_data) * 0.2), 15)
  } else {
    pdf_width <- 8
    pdf_height <- 6
  }
  
  pdf(heat_map, width = pdf_width, height = pdf_height)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  
  cat(sprintf("Heatmap saved to %s\n", heat_map))
}

# Main execution
motif_profile <- "/home/shuaiw/borg/paper/motif_change/result_drep2_99/248_1.csv"
heat_map <- "../../tmp/figures/strain_diff/248_1_heatmap.pdf"

# Define motifs to remove
motifs_to_remove <- c("GATCWADWD_2", "GATCAAGD_2", "TGATCAKWND_3", "DGATCWWDND_3", 
                      "HGATCNHHTA_3", "GATCAWDHD_2", "WGATCMNWH_3", "GATCHWW_2", 
                      "GATCWWD_2", "TNDGATCWNDND_5", "GATCWWSNW_2", "DTGGCCAWNV_5", 
                      "DTGGCCAH_5", "GNTGGCCADC_6", "KNTGGCCANCV_6", "GCYAAGRNA_6", 
                      "GCYAAGGYB_6", "GATCNNW_2", "DGATCNNNHW_3", "WGATCNNSW_3", 
                      "VAAACCGGY_6")

# Read data
df <- read.csv(motif_profile)

cat("Data shape:", nrow(df), "x", ncol(df), "\n")
cat("Columns:", paste(colnames(df), collapse=", "), "\n")
cat("Sample data:\n")
print(head(df))

# Create heatmap
heatmap_plot(df, heat_map, motifs_to_remove)
