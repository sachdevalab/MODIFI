#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)

plot_cross_taxa <- function(fig_dir) {
  # Read data
  jaccard_all <- read.csv(paste0(fig_dir, "/jaccard_all_samples.csv"))
  

  ## replace same_sample to same_isolate
  # jaccard_all$relation[jaccard_all$relation == "same_sample"] <- "same_isolate"
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
              max_y = max(jaccard_similarity_filtered2, na.rm = TRUE),
              .groups = 'drop')
  
  # Calculate mean and standard error for each relation
  stats_data <- jaccard_all %>%
    group_by(relation) %>%
    summarise(
      mean = mean(jaccard_similarity_filtered2, na.rm = TRUE),
      se = sd(jaccard_similarity_filtered2, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  # Calculate proportion with perfect match (jaccard == 1)
  perfect_match_data <- jaccard_all %>%
    group_by(relation) %>%
    summarise(
      count = n(),
      perfect_count = sum(jaccard_similarity_filtered2 == 1, na.rm = TRUE),
      proportion = perfect_count / count,
      .groups = 'drop'
    )
  
  # Calculate proportion with perfect match (jaccard == 1) for all data
  perfect_match_data <- jaccard_all %>%
    group_by(relation) %>%
    summarise(
      count = n(),
      perfect_count = sum(jaccard_similarity_filtered2 == 1, na.rm = TRUE),
      proportion = perfect_count / count,
      .groups = 'drop'
    )
  
  # calculate proportion with perfect match (jaccard_similarity_filtered == 1) for all data
  perfect_match_filtered_data <- jaccard_all %>%
    group_by(relation) %>%
    summarise(
      count = n(),
      perfect_count = sum(jaccard_similarity_filtered == 1, na.rm = TRUE),
      proportion = perfect_count / count,
      .groups = 'drop'
    )

  # Print statistics
  cat("\n=== Mean Jaccard Similarity by Relation ===\n")
  print(stats_data %>% select(relation, mean, se), n = Inf)
  
  cat("\n=== Proportion with Perfect Match (Jaccard == 1) - All Data ===\n")
  print(perfect_match_data %>% mutate(proportion_pct = round(proportion * 100, 2)), n = Inf)
  
  cat("\n=== Proportion with Perfect Match (jaccard_similarity_filtered == 1) ===\n")
  print(perfect_match_filtered_data %>% mutate(proportion_pct = round(proportion * 100, 2)), n = Inf)
  
  cat("\n=== Same Isolate with Lowest Jaccard Similarity ===\n")
  same_isolate_data <- jaccard_all %>%
    filter(relation == 'same_isolate') %>%
    arrange(jaccard_similarity_filtered2) %>%
    head(20) %>%
    select(mge_contig, host_contig, jaccard_similarity2, jaccard_similarity_filtered2, 
           all_num_filtered, share_num, mge_type, mge_length, host_length)
  print(as.data.frame(same_isolate_data), row.names = FALSE)
  cat("\n")

  # Create barplot with error bars
  p2 <- ggplot(stats_data, aes(x = relation, y = mean, fill = relation)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                  width = 0.2, linewidth = 0.5) +
    geom_text(data = count_data, 
              aes(x = relation, y = -max(stats_data$mean) * 0.05, label = paste0("n=", count)),
              size = 3.5, inherit.aes = FALSE) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      x = "MGE-Genome Taxonomic Relationship",
      y = "Mean Jaccard Similarity",
      title = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  # Save plot
  ggsave(paste0(fig_dir, "/jaccard_similarity_by_taxa_relation.pdf"), 
         p2, width = 6, height = 5.5, dpi = 400)
  
  cat("\nPlot saved successfully!\n")
}

# Main execution
fig_dir <- "../../tmp/figures/motif_sharing"

# Run analysis
plot_cross_taxa(fig_dir)
