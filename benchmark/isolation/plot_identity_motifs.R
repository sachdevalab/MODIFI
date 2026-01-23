#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)

analyze_jaccard <- function(fig_dir) {
  # Read data

  same_sample_df <- read.csv(paste0(fig_dir, "/jaccard_same_sample.csv"))
  
  # Calculate overall proportion for jaccard_similarity
  original_match <- same_sample_df %>% filter(jaccard_similarity == 1)
  proportion_original <- nrow(original_match) / nrow(same_sample_df)
  
  cat(sprintf("\nOverall proportion with perfect Jaccard similarity (=1): %.4f (%d/%d)\n",
              proportion_original, nrow(original_match), nrow(same_sample_df)))

  # Calculate proportion of rows with jaccard_similarity_filtered == 1
  perfect_match <- same_sample_df %>% filter(jaccard_similarity_filtered == 1)
  proportion_perfect <- nrow(perfect_match) / nrow(same_sample_df)
  
  cat(sprintf("Overall proportion with perfect Jaccard similarity_filtered (=1): %.4f (%d/%d)\n",
              proportion_perfect, nrow(perfect_match), nrow(same_sample_df)))


  # Calculate proportion by phylum for jaccard_similarity
  phylum_df1 <- same_sample_df %>%
    group_by(phylum) %>%
    summarise(
      perfect_count = sum(jaccard_similarity == 1),
      total_count = n(),
      proportion = perfect_count / total_count,
      .groups = 'drop'
    )
  
  cat("\n--- By Phylum (jaccard_similarity) ---\n")
  cat(sprintf("Total: %.4f (%d/%d)\n", 
              sum(phylum_df1$perfect_count) / sum(phylum_df1$total_count),
              sum(phylum_df1$perfect_count), sum(phylum_df1$total_count)))
  cat("\nProportion of perfect Jaccard similarity (=1) by phylum:\n")
  for (i in 1:nrow(phylum_df1)) {
    cat(sprintf("%s: %.4f (%d/%d)\n", 
                phylum_df1$phylum[i], 
                phylum_df1$proportion[i],
                phylum_df1$perfect_count[i],
                phylum_df1$total_count[i]))
  }
  
  # Plot bar plot for jaccard_similarity
  p1 <- ggplot(phylum_df1, aes(x = phylum, y = proportion)) +
    geom_bar(stat = "identity", fill = "gray60", width = 0.7) +
    geom_text(aes(label = paste0("n=", total_count)), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "Phylum",
      y = "Proportion of MGE-host pairs with identical motif sets",
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank()
    )
  
  # Calculate proportion by phylum for jaccard_similarity_filtered
  phylum_df2 <- same_sample_df %>%
    group_by(phylum) %>%
    summarise(
      perfect_count = sum(jaccard_similarity_filtered == 1),
      total_count = n(),
      proportion = perfect_count / total_count,
      .groups = 'drop'
    )
  
  cat("\n--- By Phylum (jaccard_similarity_filtered) ---\n")
  cat(sprintf("Total: %.4f (%d/%d)\n", 
              sum(phylum_df2$perfect_count) / sum(phylum_df2$total_count),
              sum(phylum_df2$perfect_count), sum(phylum_df2$total_count)))
  cat("\nProportion of perfect Jaccard similarity_filtered (=1) by phylum:\n")
  for (i in 1:nrow(phylum_df2)) {
    cat(sprintf("%s: %.4f (%d/%d)\n", 
                phylum_df2$phylum[i], 
                phylum_df2$proportion[i],
                phylum_df2$perfect_count[i],
                phylum_df2$total_count[i]))
  }
  
  # Plot bar plot for jaccard_similarity_filtered
  p2 <- ggplot(phylum_df2, aes(x = phylum, y = proportion)) +
    geom_bar(stat = "identity", fill = "gray60", width = 0.7) +
    geom_text(aes(label = paste0("n=", total_count)), 
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "Phylum",
      y = "Proportion with identical motifs (MGE-validated sites only)",
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank()
    )
  
  # Combine plots and save
  combined_plot <- arrangeGrob(p1, p2, ncol = 2)
  ggsave(paste0(fig_dir, "/proportion_perfect_jaccard_by_phylum.pdf"), 
         combined_plot, width = 7, height = 5, dpi = 400)
  
  cat("\nPlots saved successfully!\n")
}

# Main execution
fig_dir <- "../../tmp/figures/motif_sharing"

# Run analysis
analyze_jaccard(fig_dir)
