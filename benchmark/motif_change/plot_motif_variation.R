#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

plot_variation_fraction <- function(data_file, paper_fig_dir) {
  # Read the CSV file
  df <- read.csv(data_file)
  
  # Check if data file exists and has required columns
  if (!"cutoff" %in% colnames(df) || !"proportion" %in% colnames(df)) {
    stop("Data file must contain 'cutoff' and 'proportion' columns")
  }
  
  # Plot 1: Proportion barplot
  p1 <- ggplot(df, aes(x = factor(cutoff), y = proportion)) +
    geom_bar(stat = "identity", fill = "skyblue", width = 0.7) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "",
         x = "Strain member cutoff",
         y = "Proportion of strains\nwith motif variation") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 11),
          panel.grid.major.x = element_blank(),
          plot.margin = margin(10, 10, 10, 10))
  
  # Prepare data for stacked barplot
  df_long <- df %>%
    select(cutoff, variation_clusters, no_variation_clusters) %>%
    pivot_longer(cols = c(variation_clusters, no_variation_clusters),
                 names_to = "type",
                 values_to = "count") %>%
    mutate(type = factor(type, 
                        levels = c("variation_clusters", "no_variation_clusters"),
                        labels = c("Variation", "No Variation")))
  
  # Plot 2: Stacked barplot
  p2 <- ggplot(df_long, aes(x = factor(cutoff), y = count, fill = type)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(values = c("Variation" = "skyblue", "No Variation" = "lightgray")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "",
         x = "Strain member cutoff",
         y = "No. of strains",
         fill = "") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 11),
          panel.grid.major.x = element_blank(),
          legend.position = c(0.75, 0.9),
          legend.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(10, 10, 10, 10))
  
  # Combine plots side by side
  combined_plot <- arrangeGrob(p1, p2, ncol = 2)
  
  # Save the combined plot
  ggsave(paste0(paper_fig_dir, "/motif_variation_combined.pdf"), 
         combined_plot, 
         width = 6, 
         height = 3, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/motif_variation_combined.pdf\n", paper_fig_dir))
}

# Main execution
paper_fig_dir <- "../../tmp/figures/strain_diff/iso_drep_99"
# paper_fig_dir <- "../../tmp/figures/strain_diff/drep_99"
data_file <- paste0(paper_fig_dir, "/motif_variation_data.csv")

# Check if data file exists
if (!file.exists(data_file)) {
  stop(sprintf("Data file not found: %s\n", data_file))
}

# Run plotting
plot_variation_fraction(data_file, paper_fig_dir)
