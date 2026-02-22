#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(ggsignif)
library(readr)

count_good_ctgs <- function(df_all_data, fig_dir) {
  # Count statistics
  cat(sprintf("%d contigs in total.\n", nrow(df_all_data)))
  
  # Extract species from lineage
  df_all_data$species <- sapply(strsplit(as.character(df_all_data$lineage), ";"), function(x) tail(x, 1))
  
  # Remove 's__' prefix from species names
  df_all_data$species_clean <- ifelse(startsWith(df_all_data$species, "s__"), 
                                      substring(df_all_data$species, 4), 
                                      df_all_data$species)
  
  unique_species <- unique(df_all_data$species_clean[df_all_data$species_clean != ""])
  unique_phyla <- unique(df_all_data$phylum)
  no_species_count <- sum(df_all_data$species_clean == "")
  
  cat(sprintf("Total unique species: %d\n", length(unique_species)))
  cat(sprintf("Total unique phyla: %d\n", length(unique_phyla)))
  cat(sprintf("Total contigs with no species annotation: %d\n", no_species_count))
  
  # Count contigs with motif_num > 0
  motif_more_than_0 <- df_all_data[df_all_data$motif_num > 0, ]
  proportion <- nrow(motif_more_than_0) / nrow(df_all_data) * 100
  cat(sprintf("Total contigs with motif number > 0: %d\n", nrow(motif_more_than_0)))
  cat(sprintf("Proportion of contigs with motif number > 0: %.2f%%\n", proportion))
  
  # Print contigs with 0 motifs
  motif_0_ctgs <- df_all_data$contig[df_all_data$motif_num == 0]
  cat(sprintf("Contigs with 0 motifs: %s\n", paste(motif_0_ctgs, collapse=", ")))
  
  # Count proportions by environment
  env_summary <- df_all_data %>%
    group_by(environment) %>%
    summarise(
      total_ctgs = n(),
      motif_pos = sum(motif_num > 0),
      proportion = sum(motif_num > 0) / n() * 100
    )
  
  # Count number of unique samples per environment
  sample_counts <- df_all_data %>%
    group_by(environment) %>%
    summarise(n_samples = n_distinct(sample))
  
  cat("\n=== Sample counts by environment ===\n")
  for (i in 1:nrow(sample_counts)) {
    cat(sprintf("Environment: %s, Number of samples: %d\n",
                sample_counts$environment[i], sample_counts$n_samples[i]))
  }
  cat(sprintf("Total unique samples across all environments: %d\n\n", n_distinct(df_all_data$sample)))
  
  for (i in 1:nrow(env_summary)) {
    cat(sprintf("Environment: %s, Total contigs: %d, Contigs with motif number > 0: %d, Proportion: %.2f%%\n",
                env_summary$environment[i], env_summary$total_ctgs[i], 
                env_summary$motif_pos[i], env_summary$proportion[i]))
  }
  
  # Count proportions by phylum
  phylum_summary <- df_all_data %>%
    group_by(phylum) %>%
    summarise(
      total_ctgs = n(),
      motif_pos = sum(motif_num > 0),
      proportion = sum(motif_num > 0) / n() * 100
    )
  
  for (i in 1:nrow(phylum_summary)) {
    cat(sprintf("Phylum: %s, Total contigs: %d, Contigs with motif number > 0: %d, Proportion: %.2f%%\n",
                phylum_summary$phylum[i], phylum_summary$total_ctgs[i], 
                phylum_summary$motif_pos[i], phylum_summary$proportion[i]))
  }
  
  # Prepare data for plotting
  ctg_num_cutoff <- 40
  
  # Filter environments with >= 40 contigs
  env_proportions <- env_summary %>%
    filter(total_ctgs >= ctg_num_cutoff) %>%
    arrange(desc(proportion))
  
  # Filter phyla with >= 40 contigs and clean phylum names
  phylum_proportions <- phylum_summary %>%
    filter(total_ctgs >= ctg_num_cutoff) %>%
    mutate(phylum_clean = gsub("^p__", "", phylum)) %>%
    arrange(desc(proportion))
  
  # Plot 1: Environment proportions (bar plot)
  p1 <- ggplot(env_proportions, aes(x = reorder(environment, -proportion), y = proportion)) +
    geom_bar(stat = "identity", fill = "#808080", width = 0.7) +
    scale_x_discrete(labels = function(x) paste0(x, "\nn=", env_proportions$total_ctgs[match(x, env_proportions$environment)])) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "",
         x = NULL,
         y = "Proportion of genomes\n with motifs (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank())
  
  # Plot 2: Phylum proportions (bar plot) - use cleaned names
  p2 <- ggplot(phylum_proportions, aes(x = reorder(phylum_clean, -proportion), y = proportion)) +
    geom_bar(stat = "identity", fill = "#808080", width = 0.7) +
    scale_x_discrete(labels = function(x) paste0(x, "\nn=", phylum_proportions$total_ctgs[match(x, phylum_proportions$phylum_clean)])) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "",
         x = NULL,
         y = "Proportion of genomes\n with motifs (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank())
  
  # Filter data for boxplots
  df_filtered_env <- df_all_data %>%
    filter(environment %in% env_proportions$environment)
  
  # Calculate mean motif_num by environment for ordering (descending)
  env_order <- df_filtered_env %>%
    group_by(environment) %>%
    summarise(mean_motif = mean(motif_num)) %>%
    arrange(desc(mean_motif)) %>%
    pull(environment)
  
  df_filtered_env$environment <- factor(df_filtered_env$environment, levels = env_order)
  
  # Get counts for environment boxplot
  env_counts <- df_filtered_env %>%
    group_by(environment) %>%
    summarise(n = n()) %>%
    arrange(match(environment, env_order))

  # Create comparison list for ocean vs other environments
  ocean_comparisons <- lapply(setdiff(env_order, "ocean"), function(env) c("ocean", env))

  # Calculate and print p-values for ocean vs other environments
  cat("\n=== T-test p-values: Ocean vs Other Environments (Motif Number) ===\n")
  for (comp in ocean_comparisons) {
    ocean_data <- df_filtered_env$motif_num[df_filtered_env$environment == comp[1]]
    other_data <- df_filtered_env$motif_num[df_filtered_env$environment == comp[2]]
    
    if (length(ocean_data) > 0 && length(other_data) > 0) {
      test_result <- t.test(ocean_data, other_data)
      cat(sprintf("%s vs %s: p-value = %.4e (mean: %.2f vs %.2f)\n", 
                  comp[1], comp[2], test_result$p.value,
                  mean(ocean_data), mean(other_data)))
    }
  }
  
  # Calculate all pairwise comparisons for motif number and print only if p < 0.05
  cat("\n=== Additional Significant Pairwise Comparisons (Motif Number, p < 0.05) ===\n")
  env_list <- as.character(env_order)
  all_pairs <- combn(env_list, 2, simplify = FALSE)
  significant_found <- FALSE
  
  for (pair in all_pairs) {
    # Skip if it's an ocean comparison (already printed above)
    if ("ocean" %in% pair) next
    
    data1 <- df_filtered_env$motif_num[df_filtered_env$environment == pair[1]]
    data2 <- df_filtered_env$motif_num[df_filtered_env$environment == pair[2]]
    
    if (length(data1) > 0 && length(data2) > 0) {
      test_result <- t.test(data1, data2)
      if (test_result$p.value < 0.05) {
        cat(sprintf("%s vs %s: p-value = %.4e (mean: %.2f vs %.2f)\n", 
                    pair[1], pair[2], test_result$p.value,
                    mean(data1), mean(data2)))
        significant_found <- TRUE
      }
    }
  }
  
  if (!significant_found) {
    cat("No additional significant pairwise comparisons found (p < 0.05)\n")
  }
  cat("\n")

  # Plot 3: Environment motif number distribution (box plot)
  p3 <- ggplot(df_filtered_env, aes(x = environment, y = motif_num)) +
    geom_boxplot(fill = "#A0A0A0", outlier.size = 0.5, width = 0.6) +
    geom_signif(comparisons = ocean_comparisons,
                test = "t.test",
                map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                step_increase = 0.12,
                tip_length = 0.02,
                textsize = 3.5) +
    scale_x_discrete(labels = function(x) paste0(x, "\nn=", env_counts$n[match(x, env_counts$environment)])) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = "",
         x = NULL,
         y = "No. of motifs per genome") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank())
  
  
  # Filter data for phylum boxplot and clean phylum names
  df_filtered_phylum <- df_all_data %>%
    filter(phylum %in% phylum_proportions$phylum) %>%
    mutate(phylum_clean = gsub("^p__", "", phylum))
  
  # Calculate mean motif_num by phylum for ordering (descending)
  phylum_order <- df_filtered_phylum %>%
    group_by(phylum_clean) %>%
    summarise(mean_motif = mean(motif_num)) %>%
    arrange(desc(mean_motif)) %>%
    pull(phylum_clean)
  
  df_filtered_phylum$phylum_clean <- factor(df_filtered_phylum$phylum_clean, levels = phylum_order)
  
  # Get counts for phylum boxplot
  phylum_counts <- df_filtered_phylum %>%
    group_by(phylum_clean) %>%
    summarise(n = n()) %>%
    arrange(match(phylum_clean, phylum_order))
  
  # Plot 4: Phylum motif number distribution (box plot)
  p4 <- ggplot(df_filtered_phylum, aes(x = phylum_clean, y = motif_num)) +
    geom_boxplot(fill = "#A0A0A0", outlier.size = 0.5, width = 0.6) +
    scale_x_discrete(labels = function(x) paste0(x, "\nn=", phylum_counts$n[match(x, phylum_counts$phylum_clean)])) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(title = "",
         x = NULL,
         y = "Motif Number") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8.5),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank())
  
  # Combine plots in 2x2 layout
  combined_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  
  # Calculate proportion of genomes with motifs per sample grouped by environment
  sample_proportions <- df_all_data %>%
    group_by(environment, sample) %>%
    summarise(
      total_genomes = n(),
      genomes_with_motifs = sum(motif_num > 0),
      proportion = (genomes_with_motifs / total_genomes) * 100,
      .groups = "drop"
    ) %>%
    filter(total_genomes >= 5)
  
  # Filter environments with sufficient data
  envs_to_include <- env_summary %>%
    filter(total_ctgs >= ctg_num_cutoff) %>%
    pull(environment)
  
  sample_proportions_filtered <- sample_proportions %>%
    filter(environment %in% envs_to_include)

  ## store sample_proportions_filtered to a csv
  readr::write_csv(sample_proportions_filtered, file.path(fig_dir, "sample_proportions_filtered.csv"))
  
  # Order environments by median proportion
  env_order_box <- sample_proportions_filtered %>%
    group_by(environment) %>%
    summarise(median_prop = median(proportion)) %>%
    arrange(desc(median_prop)) %>%
    pull(environment)
  
  sample_proportions_filtered$environment <- factor(sample_proportions_filtered$environment, 
                                                     levels = env_order_box)
  
  # Create comparison list for all pairwise comparisons between environments
  # Exclude ocean as it may have only one observation
  envs_to_compare <- setdiff(as.character(env_order_box), "ocean")
  env_pairs <- combn(envs_to_compare, 2, simplify = FALSE)
  pairwise_comparisons <- list()
  
  for (pair in env_pairs) {
    n1 <- sum(sample_proportions_filtered$environment == pair[1])
    n2 <- sum(sample_proportions_filtered$environment == pair[2])
    if (n1 >= 2 && n2 >= 2) {
      pairwise_comparisons[[length(pairwise_comparisons) + 1]] <- c(pair[1], pair[2])
    }
  }
  
sample_proportions_filtered

  # Calculate and print p-values for pairwise environment comparisons
  cat("\n=== T-test p-values: Pairwise Environment Comparisons (Proportion per Sample) ===\n")
  for (pair in pairwise_comparisons) {
    data1 <- sample_proportions_filtered$proportion[sample_proportions_filtered$environment == pair[1]]
    data2 <- sample_proportions_filtered$proportion[sample_proportions_filtered$environment == pair[2]]
    
    if (length(data1) >= 2 && length(data2) >= 2) {
      test_result <- t.test(data1, data2)
      cat(sprintf("%s (n=%d) vs %s (n=%d): p-value = %.4e (mean: %.2f%% vs %.2f%%)\n", 
                  pair[1], length(data1), pair[2], length(data2),
                  test_result$p.value, mean(data1), mean(data2)))
    }
  }
  cat("\n")

  # Create boxplot of proportion per sample by environment with t.test
  p5 <- ggplot(sample_proportions_filtered, aes(x = environment, y = proportion)) +
    geom_boxplot(fill = "#A0A0A0", outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "#404040") +
    # Highlight samples with proportion == 0
    geom_point(data = subset(sample_proportions_filtered, proportion == 0),
           aes(x = environment, y = proportion),
           color = "red", size = 3, shape = 17, inherit.aes = FALSE) +
    # Label those samples
    geom_text(data = subset(sample_proportions_filtered, proportion == 0),
          aes(x = environment, y = proportion, label = sample),
          vjust = 1.5, color = "red", size = 2.5, angle = 45, inherit.aes = FALSE)
  
    # Add significance marks at the top, split lines for each comparison
    if (length(pairwise_comparisons) > 0) {
      # Stack the lines at the top, e.g. 103, 108, 113, ...
      y_positions <- seq(103, 103 + (length(pairwise_comparisons) - 1) * 10, by = 10)
      p5 <- p5 + geom_signif(comparisons = pairwise_comparisons,
                             test = "t.test",
                             map_signif_level = c("***" = 0.001, "**" = 0.01, "*" = 0.05),
                             y_position = y_positions,
                             tip_length = 0.02,
                             textsize = 3,
                             na.rm = TRUE)
    }
  
  p5 <- p5 +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0)),
        limits = c(0, 130),
        breaks = seq(0, 100, 25),
        labels = seq(0, 100, 25),
        oob = scales::squish
      ) +
      coord_cartesian(ylim = c(0, 130), clip = "off") +
    labs(title = "",
         x = NULL,
         y = "Proportion of genomes\nwith motifs per sample (%)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
          panel.grid.major.x = element_blank())
  
  # Save the combined plot
  ggsave(paste0(fig_dir, "/proportion_motif.pdf"), 
         combined_plot, 
         width = 8, 
         height = 5, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/proportion_motif.pdf\n", fig_dir))
  
  # Save p1 separately
  ggsave(paste0(fig_dir, "/proportion_motif_env_barplot.pdf"), 
         p1, 
         width = 3, 
         height = 3, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/proportion_motif_env_barplot.pdf\n", fig_dir))

  # Save p2 separately
  ggsave(paste0(fig_dir, "/proportion_motif_phylum_barplot.pdf"), 
         p2, 
         width = 6, 
         height = 3, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/proportion_motif_phylum_barplot.pdf\n", fig_dir))
  
  # Save p3 separately
  ggsave(paste0(fig_dir, "/motif_num_env_boxplot.pdf"), 
         p3, 
         width = 3, 
         height = 4, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/motif_num_env_boxplot.pdf\n", fig_dir))
  
  # Save p4 separately
  ggsave(paste0(fig_dir, "/motif_num_phylum_boxplot.pdf"), 
         p4, 
         width = 6, 
         height = 3, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/motif_num_phylum_boxplot.pdf\n", fig_dir))
  
  # Save p5 separately
  ggsave(paste0(fig_dir, "/proportion_per_sample_by_env_boxplot.pdf"), 
         p5, 
         width = 3, 
         height = 4, 
         dpi = 300)
  
  cat(sprintf("Plot saved to %s/proportion_per_sample_by_env_boxplot.pdf\n", fig_dir))
}

# Main execution
# Read data
fig_dir <- "../../tmp/figures/multi_env_linkage/"
df_all_data <- read.csv(paste0(fig_dir, "motif_num_all_samples.csv"))

# Run analysis
count_good_ctgs(df_all_data, fig_dir)

# Count the average number of motifs per genome
average_motifs <- mean(df_all_data$motif_num)
cat(sprintf("Average number of motifs per genome: %.2f, median: %.2f\n", average_motifs, median(df_all_data$motif_num)))
