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
  
  # Define labels with same_strain having two lines
  relation_labels <- c('same_isolate', 'same_strain\n(different_isolation)', 'same_species\n(different_strain)', 'same_genus\n(different_species)', 
                       'same_family\n(different_genus)', 'same_order\n(different_family)', 'same_class\n(different_order)', 'same_phylum\n(different_class)', 'same_domain\n(different_phylum)')
  
  # Convert relation to ordered factor
  jaccard_all$relation <- factor(jaccard_all$relation, levels = relation_order, labels = relation_labels)
  
  # Print total pairs
  cat(sprintf("Total pairs: %d\n", nrow(jaccard_all)))
  
  # Calculate counts for each relation
  count_data <- jaccard_all %>%
    group_by(relation) %>%
    summarise(count = n(), 
              max_y = max(jaccard_similarity2, na.rm = TRUE),
              .groups = 'drop')



  # Create boxplot
  p <- ggplot(jaccard_all, aes(x = relation, y = jaccard_similarity2, fill = relation)) +
    geom_boxplot(show.legend = FALSE) +
    geom_text(data = count_data, 
              aes(x = relation, y = max_y, label = paste0("n=", count)),
              vjust = -0.5, size = 3.5, inherit.aes = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "MGE-Genome Taxonomic Relationship",
      y = "Motif Jaccard Similarity",
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
         p, width = 6, height = 5.5, dpi = 400)
  
  cat("\nPlot saved successfully!\n")
}

# Main execution
fig_dir <- "../../tmp/figures/motif_sharing"

# Run analysis
plot_cross_taxa(fig_dir)
