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
  
  # Remove edge weights if they exist
  if ("weight" %in% edge_attr_names(G)) {
    G <- delete_edge_attr(G, "weight")
  }

  cat("Network loaded:\n")
  cat("  Nodes:", vcount(G), "\n")
  cat("  Edges:", ecount(G), "\n")

  # Get connected components
  # Use igraph::decompose to get all components as subgraphs
  comp_list <- decompose(G)
  # Only keep components with >5 nodes
  comp_list <- Filter(function(g) vcount(g) > 5, comp_list)
  # Sort by size descending
  comp_list <- comp_list[order(sapply(comp_list, vcount), decreasing = TRUE)]
  plots <- list()

  # --- GLOBAL COLOR/SHAPE MAPS ---
  all_types <- unique(V(G)$type)
  host_types <- setdiff(all_types, c("virus", "plasmid", "novel"))
  host_palette <- NULL
  if (length(host_types) > 0) {
    host_palette <- colorRampPalette(brewer.pal(12, "Set3"))(length(host_types))
    names(host_palette) <- host_types
  }
  color_map <- setNames(rep("gray", length(all_types)), all_types)
  shape_map <- setNames(rep(22, length(all_types)), all_types)
  color_map["virus"] <- "red"
  color_map["plasmid"] <- "blue"
  color_map["novel"] <- "orange"
  shape_map["virus"] <- 24
  shape_map["plasmid"] <- 21
  shape_map["novel"] <- 23
  if (!is.null(host_palette)) {
    color_map[names(host_palette)] <- host_palette
    shape_map[names(host_palette)] <- 22
  }

  # --- Circular layout for each component, arranged in a single panel ---
  # Compute circular layouts and offset each component
  layout_list <- list()
  node_df_list <- list()
  edge_df_list <- list()
  x_offset <- 0
  y_offset <- 0
  ncol_grid <- 6  # Number of columns in the grid arrangement
  grid_ncol <- ncol_grid
  grid_spacing <- 8  # Spacing between component centers
    nrow_grid <- ceiling(length(comp_list) / ncol_grid)
  for (i in seq_along(comp_list)) {
    subG <- comp_list[[i]]
    node_types <- V(subG)$type
    if ("name" %in% vertex_attr_names(subG)) {
      node_names <- V(subG)$name
    } else {
      node_names <- as.character(1:vcount(subG))
    }
    nodes_df <- data.frame(
      name = node_names,
      type = node_types,
      color = color_map[node_types],
      shape = as.numeric(shape_map[node_types]),
      stringsAsFactors = FALSE
    )
    V(subG)$color <- nodes_df$color
    V(subG)$shape <- nodes_df$shape

    # Circular layout
    circ_layout <- layout_in_circle(subG)
    # Center the layout at (x_offset, y_offset)
    circ_layout[,1] <- circ_layout[,1] + x_offset
    circ_layout[,2] <- circ_layout[,2] + y_offset
    nodes_df$x <- circ_layout[,1]
    nodes_df$y <- circ_layout[,2]
    nodes_df$comp <- i
    layout_list[[i]] <- circ_layout
    node_df_list[[i]] <- nodes_df

    # Edges for plotting
    edges <- as.data.frame(igraph::as_data_frame(subG, what = "edges"))
    if (nrow(edges) > 0) {
      edges$x <- circ_layout[match(edges$from, node_names), 1]
      edges$y <- circ_layout[match(edges$from, node_names), 2]
      edges$xend <- circ_layout[match(edges$to, node_names), 1]
      edges$yend <- circ_layout[match(edges$to, node_names), 2]
      edge_df_list[[i]] <- edges
    }

    # Update offsets for grid arrangement
    if (i %% grid_ncol == 0) {
      x_offset <- 0
      y_offset <- y_offset - grid_spacing
    } else {
      x_offset <- x_offset + grid_spacing
    }
  }

  all_nodes <- do.call(rbind, node_df_list)
  all_edges <- do.call(rbind, edge_df_list)

  # Plot all components in a single ggplot panel
  library(ggplot2)
  library(cowplot)
  p <- ggplot() +
    geom_segment(data = all_edges, aes(x = x, y = y, xend = xend, yend = yend), color = "gray", alpha = 0.7, linewidth = 0.4) +
    geom_point(data = all_nodes, aes(x = x, y = y, fill = color, shape = shape), size = 5, color = "black", stroke = 0.3, alpha = 0.95) +
    scale_fill_identity(name = "Type/Phylum", guide = "legend", labels = names(color_map), breaks = unname(color_map)) +
    scale_shape_identity(name = "Type/Phylum", guide = "legend", labels = names(shape_map), breaks = as.numeric(shape_map)) +
    theme_void() +
    theme(legend.position = "bottom", legend.box = "horizontal")

  # Legend extraction as before
  dummy_types <- names(color_map)
  dummy_df <- data.frame(
    color = unname(color_map),
    shape = as.numeric(shape_map),
    type = dummy_types
  )
  dummy_plot <- ggplot(dummy_df, aes(x = type, y = 1, fill = color, shape = shape)) +
    geom_point(size = 6, show.legend = TRUE) +
    scale_fill_identity(name = "Type/Phylum", guide = "legend", labels = dummy_types, breaks = unname(color_map)) +
    scale_shape_identity(name = "Type/Phylum", guide = "legend", labels = dummy_types, breaks = as.numeric(shape_map)) +
    theme_void() +
    theme(legend.position = "bottom", legend.box = "horizontal")
  get_legend <- function(myplot) {
    tmp <- cowplot::plot_grid(myplot, ncol = 1)
    g <- ggplotGrob(tmp)
    legend_index <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
    if (length(legend_index) == 0) return(NULL)
    g$grobs[[legend_index]]
  }
  legend <- get_legend(dummy_plot)
  if (is.null(legend)) {
    warning("No legend could be extracted. The plot will be saved without a legend.")
    legend <- NULL
  }

  # Remove legends from all plots
  plots_nolegend <- lapply(plots, function(p) p + theme(legend.position = "none"))

  # Pad plot list to fill the grid
  n_plots <- length(plots_nolegend)
  total_slots <- nrow_grid * ncol_grid
  if (n_plots < total_slots) {
    empty_plot <- function() ggplot() + theme_void()
    plots_nolegend <- c(plots_nolegend, replicate(total_slots - n_plots, empty_plot(), simplify = FALSE))
  }

  # Combine plots
  combined_plot <- plot_grid(plotlist = plots_nolegend, ncol = ncol_grid, align = "hv", rel_widths = rep(1, ncol_grid), rel_heights = rep(1, nrow_grid), hjust = -0.1, vjust = 1.1)

  # Add legend if available
  if (!is.null(legend)) {
    final_plot <- plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.08))
  } else {
    final_plot <- combined_plot
  }

  cat("Saving plot to:", output_file, "\n")
  # Set width and height for a compact, publication-ready layout
    plot_width <- 5 * ncol_grid  
    plot_height <- 5 * nrow_grid 
  plot_width <- min(plot_width, 20)
  plot_height <- min(plot_height, 20)
  ggsave(output_file, plot = final_plot, width = plot_width, height = plot_height, dpi = 600, limitsize = FALSE)
  cat("Done!\n")
}

# Main execution
if (!interactive()) {
  paper_fig_dir <- "../../tmp/figures/multi_env_linkage/network_99"
  gml_file <- file.path(paper_fig_dir, "whole_network.gml")
  output_file <- file.path(paper_fig_dir, "network_plot.pdf")
  
  plot_network(gml_file, output_file)
}
