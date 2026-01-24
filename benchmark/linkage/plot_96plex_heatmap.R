library(ComplexHeatmap)
library(dplyr)
library(readr)
library(circlize)
library(grid)
library(tidyr)
library(tibble)
library(RColorBrewer)

# Function to read and process long format data
process_long_data <- function(motif_profile_long) {
  # Read long format data
  df_long <- read_csv(motif_profile_long, show_col_types = FALSE)
  
  # Convert to wide format (pivot)
  df_wide <- df_long %>%
    select(motif_identifier, sample, fraction) %>%
    pivot_wider(names_from = sample, values_from = fraction)
  
  # Set motif_identifier as row names
  motif_mat <- df_wide %>%
    column_to_rownames("motif_identifier") %>%
    as.matrix()
  
  # Get strain and genome type info for each sample
  sample_info <- df_long %>%
    select(sample, strain_name, genome_type) %>%
    distinct() %>%
    arrange(sample)
  
  return(list(matrix = motif_mat, sample_info = sample_info))
}

# Function to create heatmap with strain annotation
create_heatmap <- function(mat, sample_info, heat_map) {
  # Remove specific motifs
  motifs_to_remove <- c("GATCBAMY_3", "BGATCGRAB_4", "GATCNAVYRNB_3", "SAAAGAGMH_6",
                         "CGWCGVKRD_2", "HCCTTCYNDD_2", "DNCCTTCYNDD_3", "GAAGGYBKND_4",
                         "CGWCGNNY_2", "RCTCGAGNRD_2", "GWCAYH_4", "BNNCCGGNYRDNA_4", 
                         "HGCGCGCNYD_3")
  mat <- mat[!rownames(mat) %in% motifs_to_remove, ]
  
  # Remove motifs with all zero values across all samples
  row_sums <- rowSums(mat)
  mat <- mat[row_sums > 0, ]
  cat("Retained", nrow(mat), "motifs with non-zero values\n")
  
  # Transpose matrix (rows=samples, cols=motifs)
  mat_t <- t(mat)
  
  # Convert values lower than 0.3 to 0
  mat_t[mat_t < 0.4] <- 0
  
  # Match sample_info to matrix rows
  sample_info <- sample_info[match(rownames(mat_t), sample_info$sample), ]
  
  # Get unique strains and assign numbers
  unique_strains <- unique(sample_info$strain_name)
  n_strains <- length(unique_strains)
  
  # Create numeric mapping for strains
  strain_numbers <- setNames(1:n_strains, unique_strains)
  sample_info$strain_number <- as.character(strain_numbers[sample_info$strain_name])
  
  # Create better color palette for strain numbers using RColorBrewer
  # Use a combination of Set3, Pastel1, and Set2 for distinct, pleasant colors
  if (n_strains <= 12) {
    strain_colors_palette <- RColorBrewer::brewer.pal(min(12, max(3, n_strains)), "Set3")
  } else if (n_strains <= 20) {
    strain_colors_palette <- c(
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(min(8, n_strains - 12), "Pastel1")
    )
  } else {
    strain_colors_palette <- c(
      RColorBrewer::brewer.pal(12, "Set3"),
      RColorBrewer::brewer.pal(9, "Pastel1"),
      RColorBrewer::brewer.pal(8, "Set2")
    )
  }
  
  strain_number_colors <- setNames(
    strain_colors_palette[1:n_strains],
    as.character(1:n_strains)
  )
  
  # Create row annotations
  row_ha <- rowAnnotation(
    Strain = sample_info$strain_number,
    Genome = sample_info$genome_type,
    col = list(
      Strain = strain_number_colors,
      Genome = c(plasmid = '#9370DB', host = '#4CAF50', unknown = '#BDBDBD')
    ),
    annotation_width = unit(c(8, 8), "mm"),
    annotation_legend_param = list(
      Genome = list(
        title = "Genome",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = unit(5, "mm"),
        grid_width = unit(5, "mm")
      ),
      Strain = list(
        title = "Strain",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels = rep("", n_strains),
        grid_height = unit(5, "mm"),
        grid_width = unit(5, "mm"),
        nrow = ceiling(n_strains / 2)
      )
    ),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "top"
  )
  
  # Define color function (white to dark red/brown)
  col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#CD5C5C", "#8B0000"))
  
  # Create heatmap with better layout
  ht <- Heatmap(
    mat_t,
    name = "Fraction",
    col = col_fun,
    left_annotation = row_ha,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = 'average',
    clustering_method_columns = 'average',
    clustering_distance_rows = 'euclidean',
    clustering_distance_columns = 'euclidean',
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 90,
    row_dend_width = unit(20, "mm"),
    column_dend_height = unit(20, "mm"),
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
      title = "Fraction",
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      grid_height = unit(5, "mm"),
      grid_width = unit(5, "mm"),
      legend_direction = "vertical",
      legend_height = unit(40, "mm")
    ),
    border = TRUE,
    row_gap = unit(0, "mm"),
    column_gap = unit(0, "mm")
  )
  
  # Save to PDF with better dimensions
  pdf(heat_map, width = 10, height = 7)
  draw(ht, padding = unit(c(5, 5, 5, 15), "mm"), 
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
}

# Main execution
motif_profile_long <- "/home/shuaiw/borg/paper/linkage/pure2/m64004_210929_143746.p100/motif_profile_long.csv"
heat_map <- "../../tmp/figures/motif_sharing/96plex_motif_heatmap.pdf"

# Process long format data
result <- process_long_data(motif_profile_long)
mat <- result$matrix
sample_info <- result$sample_info

# Create heatmap
create_heatmap(mat, sample_info, heat_map)

# Print strain mapping
cat("\nStrain Number Mapping:\n")
unique_strains <- unique(sample_info$strain_name)
strain_numbers <- setNames(1:length(unique_strains), unique_strains)
for (strain in names(strain_numbers)) {
  cat(sprintf("%2d: %s\n", strain_numbers[strain], strain))
}

cat("\nHeatmap saved to:", heat_map, "\n")
cat("Number of samples:", nrow(sample_info), "\n")
cat("Number of unique strains:", length(unique(sample_info$strain_name)), "\n")
cat("Number of motifs:", nrow(mat), "\n")