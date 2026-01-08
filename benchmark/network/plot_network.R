#!/usr/bin/env Rscript

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Function to plot network with custom shapes and colors
plot_network <- function(gml_file, output_file) {
  
  # Read GML file
  cat("Reading network from:", gml_file, "\n")
  G <- read_graph(gml_file, format = "gml")
  
  cat("Network loaded:\n")
  cat("  Nodes:", vcount(G), "\n")
  cat("  Edges:", ecount(G), "\n")
  
  # Get node types
  if (!"type" %in% vertex_attr_names(G)) {
    stop("Graph does not have 'type' vertex attribute")
  }
  
  node_types <- V(G)$type
  cat("  Node types:", unique(node_types), "\n")
  
  # Count each type
  virus_count <- sum(node_types == "virus")
  plasmid_count <- sum(node_types == "plasmid")
  novel_count <- sum(node_types == "novel")
  host_count <- sum(!node_types %in% c("virus", "plasmid", "novel"))
  
  cat("  Virus nodes:", virus_count, "\n")
  cat("  Plasmid nodes:", plasmid_count, "\n")
  cat("  Unknown nodes:", novel_count, "\n")
  cat("  Host nodes:", host_count, "\n")
  
  # Prepare node data frame
  # Check if name attribute exists, otherwise use vertex index
  if ("name" %in% vertex_attr_names(G)) {
    node_names <- V(G)$name
  } else {
    node_names <- as.character(1:vcount(G))
  }
  
  nodes_df <- data.frame(
    name = node_names,
    type = node_types,
    is_host = !node_types %in% c("virus", "plasmid", "novel"),
    stringsAsFactors = FALSE
  )
  
  # For host nodes, each unique cluster gets a unique color
  # Get unique host types (phyla)
  host_types <- unique(nodes_df$type[nodes_df$is_host])
  
  # Assign colors to host clusters
  if (length(host_types) > 0) {
    # Use colorRampPalette to generate enough colors
    if (length(host_types) <= 20) {
      host_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(host_types))
    } else {
      host_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(host_types))
    }
    names(host_palette) <- host_types
  }
  
  # Assign colors to all nodes
  nodes_df$color <- sapply(1:nrow(nodes_df), function(i) {
    type <- nodes_df$type[i]
    if (type == "virus") {
      return("red")
    } else if (type == "plasmid") {
      return("blue")
    } else if (type == "novel") {
      return("orange")
    } else {
      # Host node - use phylum color
      return(host_palette[type])
    }
  })
  
  # Assign shapes
  nodes_df$shape <- sapply(nodes_df$type, function(type) {
    if (type == "virus") {
      return(23)  # filled diamond
    } else if (type == "plasmid") {
      return(21)  # filled circle
    } else if (type == "novel") {
      return(24)  # filled triangle up
    } else {
      return(22)  # filled square for hosts
    }
  })
  
  # Create layout using graphopt (similar to neato)
  cat("\nComputing layout...\n")
  set.seed(42)
  layout <- layout_with_graphopt(G)
  
  # Create ggraph plot
  cat("Creating plot...\n")
  
  p <- ggraph(G, layout = layout) +
    # Draw edges first
    geom_edge_link(color = "gray", alpha = 0.8, width = 0.5) +
    # Draw nodes
    geom_node_point(aes(color = I(nodes_df$color), shape = I(nodes_df$shape)), 
                    size = 3, alpha = 0.8) +
    # Prepare legend manually
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  
  # Create custom legend
  legend_data <- data.frame(
    label = character(),
    color = character(),
    shape = numeric(),
    stringsAsFactors = FALSE
  )
  
  if (plasmid_count > 0) {
    legend_data <- rbind(legend_data, 
                        data.frame(label = "Plasmid (○)",
                                  color = "blue", shape = 21))
  }
  if (virus_count > 0) {
    legend_data <- rbind(legend_data, 
                        data.frame(label = "Virus (◆)",
                                  color = "red", shape = 23))
  }
  if (novel_count > 0) {
    legend_data <- rbind(legend_data, 
                        data.frame(label = "Unknown (▲)",
                                  color = "orange", shape = 24))
  }
  
  # Add host label with square symbol
  if (host_count > 0) {
    legend_data <- rbind(legend_data,
                        data.frame(label = "Host (■)",
                                  color = "gray60", shape = 22))
  }
  
  # Add top 6 host types
  if (host_count > 0) {
    host_type_counts <- nodes_df %>%
      filter(is_host) %>%
      group_by(type) %>%
      summarise(count = n(), .groups = "drop") %>%
      arrange(desc(count)) %>%
      head(6)
    
    for (i in 1:nrow(host_type_counts)) {
      type <- host_type_counts$type[i]
      count <- host_type_counts$count[i]
      legend_data <- rbind(legend_data, 
                          data.frame(label = sprintf("  %s (%d)", type, count),
                                    color = host_palette[type], shape = 22))
    }
    
    remaining_host_types <- length(host_types) - min(6, nrow(host_type_counts))
    if (remaining_host_types > 0) {
      legend_data <- rbind(legend_data, 
                          data.frame(label = sprintf("  ... +%d more host types", remaining_host_types),
                                    color = "lightgray", shape = 22))
    }
  }
  
  # Create manual legend using grid and arrange below plot
  library(grid)
  library(gridExtra)
  library(cowplot)
  
  # Create legend as a separate plot
  create_legend_plot <- function(legend_data) {
    n_items <- nrow(legend_data)
    n_cols <- 3  # 3 columns for legend
    n_rows <- ceiling(n_items / n_cols)
    
    # Create empty plot for legend
    legend_plot <- ggplot() + 
      theme_void() +
      xlim(0, n_cols) +
      ylim(0, n_rows)
    
    # Add points and text for each legend item
    for (i in 1:n_items) {
      col_idx <- (i - 1) %% n_cols
      row_idx <- n_rows - floor((i - 1) / n_cols)
      
      # Add point
      legend_plot <- legend_plot +
        annotate("point", x = col_idx + 0.1, y = row_idx - 0.5,
                shape = legend_data$shape[i], 
                color = legend_data$color[i],
                fill = legend_data$color[i],
                size = 3)
      
      # Add text
      legend_plot <- legend_plot +
        annotate("text", x = col_idx + 0.2, y = row_idx - 0.5,
                label = legend_data$label[i], 
                hjust = 0, size = 3.5)
    }
    
    return(legend_plot)
  }
  
  legend_plot <- create_legend_plot(legend_data)
  
  # Combine main plot and legend
  combined_plot <- plot_grid(
    p, legend_plot,
    ncol = 1,
    rel_heights = c(1, 0.15)
  )
  
  # Save plot
  cat("Saving plot to:", output_file, "\n")
  ggsave(output_file, plot = combined_plot, width = 9, height = 9, dpi = 600)
  
  cat("Done!\n")
}

# Main execution
if (!interactive()) {
  paper_fig_dir <- "../../tmp/figures/multi_env_linkage/network_99"
  gml_file <- file.path(paper_fig_dir, "whole_network.gml")
  output_file <- file.path(paper_fig_dir, "network_plot.png")
  
  plot_network(gml_file, output_file)
}
