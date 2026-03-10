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
  
  # Add count labels to phylum names
  phylum_df1 <- phylum_df1 %>%
    mutate(phylum_label = paste0(phylum, "\nn=", total_count))
  
  # Plot bar plot for jaccard_similarity
  p1 <- ggplot(phylum_df1, aes(x = phylum_label, y = proportion)) +
    geom_bar(stat = "identity", fill = "gray60", width = 0.7) +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Proportion of MGE-host pairs with identical motifs",
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
  
  # Add count labels to phylum names
  phylum_df2 <- phylum_df2 %>%
    mutate(phylum_label = paste0(phylum, "\nn=", total_count))
  
  # Plot bar plot for jaccard_similarity_filtered
  p2 <- ggplot(phylum_df2, aes(x = phylum_label, y = proportion)) +
    geom_bar(stat = "identity", fill = "gray60", width = 0.7) +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Proportion of MGE-host pairs with identical motifs \n(motifs not occurring in MGE are removed)",
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
  
  # Calculate proportion by phylum and mge_type for jaccard_similarity
  phylum_mge_df1 <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(
      perfect_count = sum(jaccard_similarity == 1),
      total_count = n(),
      proportion = perfect_count / total_count,
      .groups = 'drop'
    )
  
  # Complete the data to include all combinations of phylum and mge_type
  all_phylums <- unique(same_sample_df$phylum)
  all_mge_types <- unique(same_sample_df$mge_type)
  all_combinations <- expand.grid(phylum = all_phylums, mge_type = all_mge_types, stringsAsFactors = FALSE)
  phylum_mge_df1 <- all_combinations %>%
    left_join(phylum_mge_df1, by = c("phylum", "mge_type")) %>%
    mutate(
      perfect_count = ifelse(is.na(perfect_count), 0, perfect_count),
      total_count = ifelse(is.na(total_count), 0, total_count),
      proportion = ifelse(is.na(proportion), 0, proportion)
    )
  
  cat("\n--- By Phylum and MGE Type (jaccard_similarity) ---\n")
  for (i in 1:nrow(phylum_mge_df1)) {
    cat(sprintf("%s - %s: %.4f (%d/%d)\n", 
                phylum_mge_df1$phylum[i],
                phylum_mge_df1$mge_type[i],
                phylum_mge_df1$proportion[i],
                phylum_mge_df1$perfect_count[i],
                phylum_mge_df1$total_count[i]))
  }
  
  # Plot bar plot by phylum and mge_type for jaccard_similarity
  p3 <- ggplot(phylum_mge_df1, aes(x = phylum, y = proportion, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = total_count), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Proportion of identical motifs",
      fill = "MGE Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
  
  # Calculate proportion by phylum and mge_type for jaccard_similarity_filtered
  phylum_mge_df2 <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(
      perfect_count = sum(jaccard_similarity_filtered == 1),
      total_count = n(),
      proportion = perfect_count / total_count,
      .groups = 'drop'
    )
  
  # Complete the data to include all combinations of phylum and mge_type
  phylum_mge_df2 <- all_combinations %>%
    left_join(phylum_mge_df2, by = c("phylum", "mge_type")) %>%
    mutate(
      perfect_count = ifelse(is.na(perfect_count), 0, perfect_count),
      total_count = ifelse(is.na(total_count), 0, total_count),
      proportion = ifelse(is.na(proportion), 0, proportion)
    )
  
  cat("\n--- By Phylum and MGE Type (jaccard_similarity_filtered) ---\n")
  for (i in 1:nrow(phylum_mge_df2)) {
    cat(sprintf("%s - %s: %.4f (%d/%d)\n", 
                phylum_mge_df2$phylum[i],
                phylum_mge_df2$mge_type[i],
                phylum_mge_df2$proportion[i],
                phylum_mge_df2$perfect_count[i],
                phylum_mge_df2$total_count[i]))
  }
  
  # Plot bar plot by phylum and mge_type for jaccard_similarity_filtered
  p4 <- ggplot(phylum_mge_df2, aes(x = phylum, y = proportion, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = total_count), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Proportion of identical motifs (filtered)",
      fill = "MGE Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
  
  # Save p1 and p2 combined
  combined_plot_12 <- arrangeGrob(p1, p2, ncol = 2)
  ggsave(paste0(fig_dir, "/proportion_perfect_jaccard_by_phylum.pdf"), 
         combined_plot_12, width = 14, height = 6, dpi = 400)
  
  # Extract legend from p3
  library(gridExtra)
  library(grid)
  get_legend <- function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  shared_legend <- get_legend(p3)
  
  # Create versions without legends
  p3_no_legend <- p3 + theme(legend.position = "none")
  p4_no_legend <- p4 + theme(legend.position = "none")
  
  # Combine p3 and p4 with shared legend at bottom
  combined_plot_34 <- arrangeGrob(
    arrangeGrob(p3_no_legend, p4_no_legend, ncol = 2),
    shared_legend,
    nrow = 2,
    heights = c(10, 0.5)
  )
  
  ggsave(paste0(fig_dir, "/proportion_perfect_jaccard_by_phylum_mge.pdf"), 
         combined_plot_34, width = 7, height = 5, dpi = 400)
  
  # Calculate mean and SE for jaccard_similarity by phylum and mge_type
  jaccard_stats <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(
      mean = mean(jaccard_similarity, na.rm = TRUE),
      se = sd(jaccard_similarity, na.rm = TRUE) / sqrt(n()),
      count = n(),
      .groups = 'drop'
    )
  
  # Complete the data to include all combinations
  jaccard_stats <- all_combinations %>%
    left_join(jaccard_stats, by = c("phylum", "mge_type")) %>%
    mutate(
      mean = ifelse(is.na(mean), 0, mean),
      se = ifelse(is.na(se), 0, se),
      count = ifelse(is.na(count), 0, count)
    )
  
    # Print mean Jaccard similarity (not filtered) by phylum and MGE type
  cat("\n--- Mean Jaccard Similarity (not filtered) by Phylum and MGE Type ---\n")
  for (i in 1:nrow(jaccard_stats)) {
    cat(sprintf("%s - %s: mean = %.4f ± %.4f (n=%d)\n", 
                jaccard_stats$phylum[i],
                jaccard_stats$mge_type[i],
                jaccard_stats$mean[i],    
                jaccard_stats$se[i],
                jaccard_stats$count[i]))
  }
  cat("\n")

  # Create barplot for jaccard_similarity by phylum and mge_type
  p5 <- ggplot(jaccard_stats, aes(x = phylum, y = mean, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.7),
                  width = 0.2, linewidth = 0.5) +
    geom_text(aes(label = count, y = mean + se), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Mean Jaccard Similarity",
      fill = "MGE Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
  
  # Calculate mean and SE for jaccard_similarity_filtered by phylum and mge_type
  jaccard_filtered_stats <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(
      mean = mean(jaccard_similarity_filtered, na.rm = TRUE),
      se = sd(jaccard_similarity_filtered, na.rm = TRUE) / sqrt(n()),
      count = n(),
      .groups = 'drop'
    )
  
  # Complete the data to include all combinations
  jaccard_filtered_stats <- all_combinations %>%
    left_join(jaccard_filtered_stats, by = c("phylum", "mge_type")) %>%
    mutate(
      mean = ifelse(is.na(mean), 0, mean),
      se = ifelse(is.na(se), 0, se),
      count = ifelse(is.na(count), 0, count)
    )
  
  # Print mean Jaccard similarity (filtered) by phylum and MGE type
  cat("\n--- Mean Jaccard Similarity (filtered) by Phylum and MGE Type ---\n")
  for (i in 1:nrow(jaccard_filtered_stats)) {
    cat(sprintf("%s - %s: mean = %.4f ± %.4f (n=%d)\n", 
                jaccard_filtered_stats$phylum[i],
                jaccard_filtered_stats$mge_type[i],
                jaccard_filtered_stats$mean[i],
                jaccard_filtered_stats$se[i],
                jaccard_filtered_stats$count[i]))
  }
  cat("\n")
  
  # Create barplot for jaccard_similarity_filtered by phylum and mge_type
  p6 <- ggplot(jaccard_filtered_stats, aes(x = phylum, y = mean, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.7),
                  width = 0.2, linewidth = 0.5) +
    geom_text(aes(label = count, y = mean + se), 
              position = position_dodge(width = 0.7), 
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(
      x = "Phylum",
      y = "Mean Jaccard Similarity (filtered)",
      fill = "MGE Type"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12),
      panel.grid.major.y = element_line(colour = "grey80"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )
  
  # Extract legend from p5 and combine p5 and p6
  shared_legend_box <- get_legend(p5)
  p5_no_legend <- p5 + theme(legend.position = "none")
  p6_no_legend <- p6 + theme(legend.position = "none")
  
  combined_plot_56 <- arrangeGrob(
    arrangeGrob(p5_no_legend, p6_no_legend, ncol = 2),
    shared_legend_box,
    nrow = 2,
    heights = c(10, 0.5)
  )
  
  ggsave(paste0(fig_dir, "/jaccard_boxplot_by_phylum_mge.pdf"), 
         combined_plot_56, width = 8.5, height = 5, dpi = 400)
  
  cat("\nPlots saved successfully!\n")
}

# Main execution
fig_dir <- "../../tmp/figures/motif_sharing"

# Run analysis
analyze_jaccard(fig_dir)
