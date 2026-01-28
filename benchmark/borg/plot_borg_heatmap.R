library(tidyr)
library(dplyr)
library(ggplot2)
library(ggdendro)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

personal_plot <- function(profile_df, cluster_species, plot_name) {
  
  # Create pivot table for heatmap
  pivot_df <- profile_df %>%
    select(contig, motifString, fraction) %>%
    pivot_wider(names_from = motifString, values_from = fraction, values_fill = 0)
  
  # Save contig names and BORG_Ref mapping BEFORE clustering
  contig_names <- pivot_df$contig
  contig_to_borg <- profile_df %>%
    group_by(contig) %>%
    summarise(BORG_Ref = first(BORG_Ref), .groups = 'drop') %>%
    tibble::deframe()
  
  # Convert to matrix for clustering (remove contig column)
  pivot_matrix <- as.matrix(pivot_df %>% select(-contig))
  rownames(pivot_matrix) <- contig_names
  
  # Create labels with BORG reference information BEFORE clustering
  # These labels match the original order of pivot_matrix rows
  y_labels_all <- sapply(contig_names, function(contig) {
    borg_ref <- contig_to_borg[contig]
    if (is.na(borg_ref)) {
      warning(paste0("Warning: ", contig, " not found in contig_to_borg mapping."))
      borg_ref <- "NA"
    }
    paste0(contig, ",", borg_ref)
  })
  
  # Extract sample names from contig names using first 8 underscore-separated parts
  extract_sample_name <- function(contig_name) {
    parts <- strsplit(contig_name, "_")[[1]]
    if (length(parts) >= 8) {
      return(paste(parts[1:8], collapse = "_"))
    } else {
      return(contig_name)
    }
  }
  
  sample_names <- sapply(contig_names, extract_sample_name)
  unique_samples <- unique(sample_names)
  
  # Assign colors to each unique sample
  n_samples <- length(unique_samples)
  # Use Dark2 palette which has better visibility (no yellow)
  if (n_samples <= 8) {
    sample_colors_palette <- brewer.pal(min(max(3, n_samples), 8), "Dark2")
  } else {
    # Combine multiple palettes for more colors, avoiding yellow
    base_colors <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))
    # Remove yellow-ish colors (index 6 in Set1 is yellow)
    base_colors <- base_colors[base_colors != "#FFFF33"]
    if (n_samples <= length(base_colors)) {
      sample_colors_palette <- base_colors[1:n_samples]
    } else {
      sample_colors_palette <- colorRampPalette(base_colors)(n_samples)
    }
  }
  names(sample_colors_palette) <- unique_samples
  
  # Create a vector of colors for each label based on sample
  label_colors <- sapply(sample_names, function(s) sample_colors_palette[s])
  
  # Perform hierarchical clustering on rows (contigs) and columns (motifs)
  # Using correlation distance and ward.D2 method (equivalent to ward in scipy)
  
  # Cluster rows (contigs)
  row_dist <- as.dist(1 - cor(t(pivot_matrix), method = "pearson"))
  row_hclust <- hclust(row_dist, method = "ward.D2")
  row_order <- row_hclust$order
  
  # Cluster columns (motifs)
  col_dist <- as.dist(1 - cor(pivot_matrix, method = "pearson"))
  col_hclust <- hclust(col_dist, method = "ward.D2")
  col_order <- col_hclust$order
  
  # Reorder pivot matrix based on clustering
  pivot_matrix_ordered <- pivot_matrix[row_order, col_order]
  
  # Create labels for tree plot (in reordered sequence)
  y_labels <- y_labels_all[row_order]
  
  # Update row names with labels
  rownames(pivot_matrix_ordered) <- y_labels
  
  # Create combined plot with dendrogram and heatmap using ComplexHeatmap
  # Define color palette for heatmap
  col_fun <- colorRamp2(seq(0, 1, length.out = 100), 
                        colorRampPalette(c("#FFFFD4", "#C7E9B4", "#41B6C4", "#225EA8", "#081D58"))(100))
  
  # Create heatmap with colored row labels
  ht <- Heatmap(pivot_matrix,
                name = "Fraction",
                col = col_fun,
                cluster_rows = row_hclust,
                cluster_columns = col_hclust,
                show_row_names = TRUE,
                show_column_names = FALSE,
                show_heatmap_legend = FALSE,
                row_labels = y_labels_all,
                row_names_gp = gpar(fontsize = 6, col = label_colors),
                column_names_gp = gpar(fontsize = 6),
                row_dend_width = unit(1.5, "cm"),
                column_dend_height = unit(1, "cm"),
                width = unit(8, "cm"))
  
  pdf(plot_name, width = 15, height = 8)
  draw(ht, padding = unit(c(2, 25, 2, 0), "mm"))
  dev.off()
  cat(paste0("Saved combined plot to: ", plot_name, "\n"))
  
  return(invisible(NULL))
}

# Example usage (uncomment and modify as needed):
profile_df <- read.csv("/home/shuaiw/borg/paper/borg_data/profile//profile_profile_df_filtered.csv")
personal_plot(profile_df, cluster_species = "borg", plot_name = "/home/shuaiw/borg/paper/borg_data/profile/borg_profile_heatmap.pdf")
