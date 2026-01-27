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
cor_test_pearson <- cor.test(dna_motif_corr_df$dnadiff, dna_motif_corr_df$jaccard, method = "pearson")
cat("Pearson P-value:", format.pval(cor_test_pearson$p.value, digits = 4), "\n")
cat("Pearson 95% Confidence Interval: [", round(cor_test_pearson$conf.int[1], 4), ",", 
    round(cor_test_pearson$conf.int[2], 4), "]\n\n")

# Format p-value for display
if (cor_test_pearson$p.value < 0.001) {
  p_text <- "P < 0.001"
} else {
  p_text <- sprintf("P = %.3f", cor_test_pearson$p.value)
}

# Create scatter plot with Pearson correlation line
p <- ggplot(dna_motif_corr_df, aes(x = dnadiff, y = jaccard)) +
  geom_point(alpha = 0.3, size = 1, color = "gray40") +
  geom_smooth(method = "lm", color = "blue", se = TRUE, alpha = 0.2, linewidth = 1) +
  annotate("text", x = max(dna_motif_corr_df$dnadiff) * 0.65, y = max(dna_motif_corr_df$jaccard) * 0.95,
           label = sprintf("Pearson r = %.3f\n%s", pearson_corr, p_text),
           hjust = 0, size = 5, fontface = "bold") +
  labs(
    x = "DNA Edit Distance",
    y = "Motif Jaccard Index"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    panel.grid.minor = element_blank()
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
  mutate(bin_label_with_n = bin_label)

# Filter bins with at least 100 samples
bin_counts_filtered <- bin_counts %>%
  filter(n >= 100)

# Add count labels to data and filter
dna_motif_corr_df <- dna_motif_corr_df %>%
  left_join(bin_counts, by = c("bin_label", "dnadiff_bin")) %>%
  filter(n >= 100)

cat("Bins after filtering (n >= 100):", nrow(bin_counts_filtered), "\n")

# Create plot with mean points and quartile error bars
# Calculate mean and quartiles for each bin
mean_jaccard_for_line <- dna_motif_corr_df %>%
  group_by(dnadiff_bin, bin_label_with_n) %>%
  summarise(
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    q25 = quantile(jaccard, 0.25, na.rm = TRUE),
    q75 = quantile(jaccard, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dnadiff_bin) %>%
  mutate(x_pos = row_number())

p_violin <- ggplot(mean_jaccard_for_line, 
                   aes(x = x_pos, y = mean_jaccard)) +
#   geom_smooth(method = "loess", color = "black", fill = "gray80", 
#               se = TRUE, linewidth = 1.2, span = 0.75, alpha = 0.3) +
  geom_line(color = "black", linewidth = 0.8, alpha = 0.5) +
  geom_errorbar(aes(ymin = q25, ymax = q75),
                width = 0.4, color = "black", linewidth = 0.8, alpha = 0.7) +
  geom_point(color = "black", size = 1.5, shape = 19) +
  annotate("text", x = max(mean_jaccard_for_line$x_pos) * 0.55, y = 0.92,
           label = sprintf("Pearson r = %.3f\n%s", pearson_corr, p_text),
           hjust = 0, size = 5, color = "black") +
  scale_x_continuous(
    breaks = mean_jaccard_for_line$x_pos,
    labels = mean_jaccard_for_line$bin_label_with_n
  ) +
  labs(
    x = "Within-strain edit distance",
    y = "Jaccard Similarity"
  ) +
  ylim(0, 1) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 13, color = "black"),
    axis.text = element_text(size = 11, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save plot
ggsave(violin_output_file, plot = p_violin, width = 4, height = 4)
cat("Saved plot to:", violin_output_file, "\n")