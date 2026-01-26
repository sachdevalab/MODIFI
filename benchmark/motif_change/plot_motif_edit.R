#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(dplyr)
library(gridExtra)

# Hard-coded paths
dna_motif_corr_file <- "../../tmp/figures/strain_diff/iso_drep_99/dna_motif_corr.csv"
paper_fig_dir <- "../../tmp/figures/strain_diff/iso_drep_99/"

output_file <- file.path(paper_fig_dir, "dna_motif_corr.pdf")
violin_output_file <- file.path(paper_fig_dir, "dna_motif_corr_violin.pdf")
combined_violin_file <- file.path(paper_fig_dir, "dna_motif_corr_violin_combined.pdf")

# Read data
dna_motif_corr_df <- read_csv(dna_motif_corr_file, show_col_types = FALSE)
  
# Check if data is empty
if (nrow(dna_motif_corr_df) == 0) {
  cat("[⚠️] No data to plot for", output_file, "\n")
  quit(status = 1)
}

# Function to find rows with edit distance = 0 but jaccard < 1
find_zero_edit_diff_jaccard <- function(df) {
  anomalies <- df %>%
    filter(dnadiff == 0 & jaccard < 1)
  
  cat("\n=== Anomaly Detection: Edit Distance = 0 but Jaccard < 1 ===\n")
  cat("Total rows in dataset:", nrow(df), "\n")
  cat("Rows with dnadiff = 0:", sum(df$dnadiff == 0), "\n")
  cat("Rows with dnadiff = 0 AND jaccard < 1:", nrow(anomalies), "\n\n")
  
  if (nrow(anomalies) > 0) {
    cat("Anomalous cases (identical DNA but different methylation):\n")
    print(anomalies)
    
    # Save to CSV
    anomaly_file <- file.path(paper_fig_dir, "anomalies_zero_edit_diff.csv")
    write_csv(anomalies, anomaly_file)
    cat("\nSaved anomalies to:", anomaly_file, "\n")
  } else {
    cat("No anomalies found.\n")
  }
  
  return(anomalies)
}

# Run anomaly detection
anomalies <- find_zero_edit_diff_jaccard(dna_motif_corr_df)

# Calculate Pearson correlation between edit distance and jaccard
cat("\n=== Correlation Analysis ===\n")
pearson_corr <- cor(dna_motif_corr_df$dnadiff, dna_motif_corr_df$jaccard, method = "pearson")
cat("Pearson correlation between DNA edit distance and Motif Jaccard Index:", 
    round(pearson_corr, 4), "\n")

# Perform correlation test for p-value
cor_test <- cor.test(dna_motif_corr_df$dnadiff, dna_motif_corr_df$jaccard, method = "pearson")
cat("P-value:", format.pval(cor_test$p.value, digits = 4), "\n")
cat("95% Confidence Interval: [", round(cor_test$conf.int[1], 4), ",", 
    round(cor_test$conf.int[2], 4), "]\n\n")

# Create scatter plot
p <- ggplot(dna_motif_corr_df, aes(x = dnadiff, y = jaccard)) +
  geom_point(alpha = 0.5) +
  labs(
    x = "DNA Edit Distance",
    y = "Motif Jaccard Index",
    title = "DNA vs Motif Correlation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Save plot
ggsave(output_file, plot = p, width = 6, height = 6)
cat("Saved plot to:", output_file, "\n")

# Bin edit distance by 5000 and create violin plot
dna_motif_corr_df <- dna_motif_corr_df %>%
  mutate(
    dnadiff_bin = floor(dnadiff / 5000) * 5000,
    bin_label = paste0(dnadiff_bin / 1000, "-", (dnadiff_bin + 5000) / 1000, "k")
  )

# Count samples in each bin
bin_counts <- dna_motif_corr_df %>%
  group_by(bin_label, dnadiff_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(dnadiff_bin) %>%
  mutate(bin_label_with_n = paste0(bin_label, "\nn=", n))

# Filter bins with at least 100 samples
bin_counts_filtered <- bin_counts %>%
  filter(n >= 100)

# Add count labels to data and filter
dna_motif_corr_df <- dna_motif_corr_df %>%
  left_join(bin_counts, by = c("bin_label", "dnadiff_bin")) %>%
  filter(n >= 100)

cat("Bins after filtering (n >= 100):", nrow(bin_counts_filtered), "\n")

# Create violin plot
p_violin <- ggplot(dna_motif_corr_df, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = jaccard)) +
  geom_violin(fill = "lightblue", alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
  labs(
    x = "DNA Edit Distance Bin (5k intervals, n≥100)",
    y = "Motif Jaccard Index",
    title = "Motif Jaccard Distribution by DNA Edit Distance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save violin plot
ggsave(violin_output_file, plot = p_violin, width = 6, height = 4)
cat("Saved violin plot to:", violin_output_file, "\n")

# Create bar plot showing mean jaccard for main bins (filtered)
mean_jaccard_main <- dna_motif_corr_df %>%
  group_by(bin_label_with_n, dnadiff_bin) %>%
  summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop")

p_bar_main <- ggplot(mean_jaccard_main, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = mean_jaccard)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  labs(
    x = "DNA Edit Distance Bin (5k intervals, n≥100)",
    y = "Mean Motif Jaccard Index",
    title = "Mean Motif Jaccard by DNA Edit Distance"
  ) +
  ylim(0, 1) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Create subplot for edit distance < 10000 with bin size 1000
dna_motif_corr_df_subset <- read_csv(dna_motif_corr_file, show_col_types = FALSE) %>%
  filter(dnadiff < 10000)

if (nrow(dna_motif_corr_df_subset) > 0) {
  # Bin edit distance by 1000
  dna_motif_corr_df_subset <- dna_motif_corr_df_subset %>%
    mutate(
      dnadiff_bin = floor(dnadiff / 1000) * 1000,
      bin_label = paste0(dnadiff_bin / 1000, "-", (dnadiff_bin + 1000) / 1000, "k")
    )
  
  # Count samples in each bin
  bin_counts_subset <- dna_motif_corr_df_subset %>%
    group_by(bin_label, dnadiff_bin) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(dnadiff_bin) %>%
    mutate(bin_label_with_n = paste0(bin_label, "\nn=", n))
  
  # Add count labels to data
  dna_motif_corr_df_subset <- dna_motif_corr_df_subset %>%
    left_join(bin_counts_subset, by = c("bin_label", "dnadiff_bin"))
  
  # Create violin plot for subset
  p_violin_subset <- ggplot(dna_motif_corr_df_subset, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = jaccard)) +
    geom_violin(fill = "lightcoral", alpha = 0.7) +
    geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
    labs(
      x = "DNA Edit Distance Bin (1k intervals)",
      y = "Motif Jaccard Index",
      title = "Motif Jaccard Distribution (Edit Distance < 10k)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Create bar plot for subset (< 10k)
  mean_jaccard_subset <- dna_motif_corr_df_subset %>%
    group_by(bin_label_with_n, dnadiff_bin) %>%
    summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop")
  
  p_bar_subset <- ggplot(mean_jaccard_subset, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = mean_jaccard)) +
    geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
    labs(
      x = "DNA Edit Distance Bin (1k intervals)",
      y = "Mean Motif Jaccard Index",
      title = "Mean Motif Jaccard (Edit Distance < 10k)"
    ) +
    ylim(0, 1) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Combine both violin plots
  combined_plot <- grid.arrange(p_violin, p_violin_subset, ncol = 1)
  ggsave(combined_violin_file, plot = combined_plot, width = 8, height = 10)
  cat("Saved combined violin plot to:", combined_violin_file, "\n")
  
  # Create another subplot for edit distance < 1000 with bin size 100
  dna_motif_corr_df_subset2 <- read_csv(dna_motif_corr_file, show_col_types = FALSE) %>%
    filter(dnadiff < 1000)
  
  if (nrow(dna_motif_corr_df_subset2) > 0) {
    # Bin edit distance by 100
    dna_motif_corr_df_subset2 <- dna_motif_corr_df_subset2 %>%
      mutate(
        dnadiff_bin = floor(dnadiff / 100) * 100,
        bin_label = paste0(dnadiff_bin, "-", dnadiff_bin + 100)
      )
    
    # Count samples in each bin
    bin_counts_subset2 <- dna_motif_corr_df_subset2 %>%
      group_by(bin_label, dnadiff_bin) %>%
      summarise(n = n(), .groups = "drop") %>%
      arrange(dnadiff_bin) %>%
      mutate(bin_label_with_n = paste0(bin_label, "\nn=", n))
    
    # Add count labels to data
    dna_motif_corr_df_subset2 <- dna_motif_corr_df_subset2 %>%
      left_join(bin_counts_subset2, by = c("bin_label", "dnadiff_bin"))
    
    # Create violin plot for subset2
    p_violin_subset2 <- ggplot(dna_motif_corr_df_subset2, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = jaccard)) +
      geom_violin(fill = "lightgreen", alpha = 0.7) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
      labs(
        x = "DNA Edit Distance Bin (100 bp intervals)",
        y = "Motif Jaccard Index",
        title = "Motif Jaccard Distribution (Edit Distance < 1k)"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Create bar plot for subset2 (< 1k)
    mean_jaccard_subset2 <- dna_motif_corr_df_subset2 %>%
      group_by(bin_label_with_n, dnadiff_bin) %>%
      summarise(mean_jaccard = mean(jaccard, na.rm = TRUE), .groups = "drop")
    
    p_bar_subset2 <- ggplot(mean_jaccard_subset2, aes(x = reorder(bin_label_with_n, dnadiff_bin), y = mean_jaccard)) +
      geom_bar(stat = "identity", fill = "lightgreen", alpha = 0.7) +
      labs(
        x = "DNA Edit Distance Bin (100 bp intervals)",
        y = "Mean Motif Jaccard Index",
        title = "Mean Motif Jaccard (Edit Distance < 1k)"
      ) +
      ylim(0, 1) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    # Combine all three violin plots
    combined_plot_all <- grid.arrange(p_violin, p_violin_subset, p_violin_subset2, ncol = 1)
    combined_all_file <- file.path(paper_fig_dir, "dna_motif_corr_violin_all.pdf")
    ggsave(combined_all_file, plot = combined_plot_all, width = 8, height = 14)
    cat("Saved all three violin plots to:", combined_all_file, "\n")
    
    # Combine all three bar plots
    combined_bar_all <- grid.arrange(p_bar_main, p_bar_subset, p_bar_subset2, ncol = 1)
    combined_bar_file <- file.path(paper_fig_dir, "dna_motif_corr_barplot_all.pdf")
    ggsave(combined_bar_file, plot = combined_bar_all, width = 8, height = 14)
    cat("Saved all three bar plots to:", combined_bar_file, "\n")
    
    # Combine violin and bar plots side by side for each scale
    combined_violin_bar <- grid.arrange(
      p_violin, p_bar_main,
      p_violin_subset, p_bar_subset,
      p_violin_subset2, p_bar_subset2,
      ncol = 2
    )
    combined_violin_bar_file <- file.path(paper_fig_dir, "dna_motif_corr_combined_all.pdf")
    ggsave(combined_violin_bar_file, plot = combined_violin_bar, width = 14, height = 16)
    cat("Saved combined violin and bar plots to:", combined_violin_bar_file, "\n")
    
  } else {
    cat("[⚠️] No data with edit distance < 1000 for third subplot\n")
  }
  
} else {
  cat("[⚠️] No data with edit distance < 10000 for subset plot\n")
}