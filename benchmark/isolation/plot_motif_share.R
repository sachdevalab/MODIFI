#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)

plot_cross_taxa <- function(fig_dir) {
  # Read data
  jaccard_all <- read.csv(paste0(fig_dir, "/jaccard_all_samples.csv"))
  

  ## replace same_sample to same_isolate
  jaccard_all$relation[jaccard_all$relation == "same_sample"] <- "same_isolate"
  # Define relation order
  relation_order <- c('same_isolate', 'same_strain', 'same_species', 'same_genus', 
                      'same_family', 'same_order', 'same_class', 'same_phylum', 'same_domain')
  
  # Convert relation to ordered factor
  jaccard_all$relation <- factor(jaccard_all$relation, levels = relation_order)
  
  # Print total pairs
  cat(sprintf("Total pairs: %d\n", nrow(jaccard_all)))
  
  # Create boxplot
  p <- ggplot(jaccard_all, aes(x = relation, y = jaccard_similarity, fill = relation)) +
    geom_boxplot(show.legend = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "MGE-Genome Taxonomic Relationship",
      y = "Modification Motif Jaccard Similarity",
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Save plot
  ggsave(paste0(fig_dir, "/jaccard_similarity_by_taxa_relation.png"), 
         p, width = 6, height = 4, dpi = 400)
  
  cat("\nPlot saved successfully!\n")
}

# Main execution
fig_dir <- "../../tmp/figures/motif_sharing"

# Run analysis
plot_cross_taxa(fig_dir)
