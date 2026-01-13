library(ComplexHeatmap)
library(dplyr)
library(readr)
library(circlize)
library(grid)

# Function to read plasmid list and get host associations
get_new_host <- function(plasmid_list_file) {
  df <- read_csv(plasmid_list_file, show_col_types = FALSE)
  
  plasmid_list <- df$seq_name
  
  # Extract and flatten host lists
  host_list <- unique(unlist(strsplit(df$host, ";")))
  
  return(list(plasmid_list = plasmid_list, host_list = host_list))
}

# Function to create heatmap
create_heatmap <- function(df, heat_map, plasmid_list, host_list) {
  # Remove specific motifs
  motifs_to_remove <- c("GATCBAMY", "BGATCGRAB", "GATCNAVYRNB")
  df <- df[!rownames(df) %in% motifs_to_remove, ]
  
  # Transpose dataframe (rows become columns, columns become rows)
  df_t <- as.data.frame(t(df))
  
  # Sort sequences
  df_t <- df_t[order(rownames(df_t)), ]
  
  # Convert to matrix for ComplexHeatmap
  mat <- as.matrix(df_t)
  
  # Convert values lower than 0.3 to 0
  mat[mat < 0.3] <- 0
  
  # Assign type to each sequence
  type_labels <- sapply(rownames(mat), function(seq) {
    if (seq %in% plasmid_list) {
      return('plasmid')
    } else if (seq %in% host_list) {
      return('host')
    } else {
      return('host')
    }
  })
  
  # Create annotation with better styling
  row_ha <- rowAnnotation(
    Type = type_labels,
    col = list(Type = c(plasmid = '#6A5ACD', host = '#FF8C00')),
    annotation_width = unit(8, "mm"),
    annotation_legend_param = list(
      Type = list(
        title = "Type", 
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 10),
        grid_height = unit(5, "mm"),
        grid_width = unit(5, "mm")
      )
    ),
    show_annotation_name = FALSE
  )
  
  # Define color function (white to dark red/brown)
  col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "#CD5C5C", "#8B0000"))
  
  # Create heatmap with better layout
  ht <- Heatmap(
    mat,
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
    row_names_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 90,
    row_dend_width = unit(20, "mm"),
    column_dend_height = unit(20, "mm"),
    rect_gp = gpar(col = "white", lwd = 1),
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
  pdf(heat_map, width = 20, height = 12)
  draw(ht, padding = unit(c(5, 5, 5, 15), "mm"), 
       heatmap_legend_side = "right",
       annotation_legend_side = "right")
  dev.off()
}

# Main execution
motif_profile <- "/home/shuaiw/borg/paper/linkage/pure2/m64004_210929_143746.p100/motif_profile.csv"
plasmid_list_file <- "/home/shuaiw/methylation/data/ZymoTrumatrix/2021-11-Microbial-96plex/ref/merged2.fa.fai.plasmid.list"
heat_map <- "../../tmp/figures/motif_sharing/96plex_motif_heatmap.pdf"

# Read motif profile - handle duplicate row names
df <- read.csv(motif_profile, check.names = FALSE)

# Check if first column should be row names
if (ncol(df) > 0 && !is.numeric(df[,1])) {
  # Make row names unique by adding suffixes to duplicates
  row_names <- make.unique(as.character(df[,1]))
  df <- df[,-1]  # Remove first column
  rownames(df) <- row_names
}

# Get plasmid and host lists
result <- get_new_host(plasmid_list_file)
plasmid_list <- result$plasmid_list
host_list <- result$host_list

# Create heatmap
create_heatmap(df, heat_map, plasmid_list, host_list)

cat("Heatmap saved to:", heat_map, "\n")