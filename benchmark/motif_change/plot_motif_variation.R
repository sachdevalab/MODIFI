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
        geom_bar(stat = "identity", fill = "gray60", width = 0.7) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(title = "",
        x = "Cluster member cutoff",
        y = "Proportion of clusters\nwith motif variation") +
      theme_minimal() +
      theme(axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 10),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(10, 10, 10, 10))

    # Save p1 as a separate file
    ggsave(paste0(paper_fig_dir, "motif_variation_proportion_gray.pdf"), p1, width = 3, height = 3, dpi = 300)
}

# Main execution
# paper_fig_dir <- "../../tmp/figures/strain_diff/iso_drep_99/"
paper_fig_dir <- "../../tmp/figures/strain_diff/drep_99/"
data_file <- paste0(paper_fig_dir, "/motif_variation_data.csv")

plot_variation_fraction(data_file, paper_fig_dir) 
