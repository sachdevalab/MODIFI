#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)
# Read data
data_file <- "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_similarity_96plex.csv"
df <- read_csv(data_file, show_col_types = FALSE)

# Keep only filtered Jaccard similarity values
df_plot <- df %>%
  select(MGE, host, jaccard_similarity_filtered, total_motif_num_filtered) %>%
  rename(
    similarity_value = jaccard_similarity_filtered,
    motif_count = total_motif_num_filtered
  )

# Reorder MGE to put S_enterica_LT2_2 first
mge_order <- unique(df_plot$MGE)
mge_order <- c("S_enterica_LT2_2", setdiff(mge_order, "S_enterica_LT2_2"))
df_plot$MGE <- factor(df_plot$MGE, levels = mge_order)

# Create the plot
p <- ggplot(df_plot, aes(x = MGE, y = similarity_value)) +
  geom_bar(stat = "identity", fill = "#5F9EA0", width = 0.7) +
  geom_text(aes(label = motif_count),
            vjust = 1.5, size = 3, color = "white") +
  labs(x = "",
       y = "Jaccard Similarity",
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 12, ),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "none") +
  ylim(0, 1)

# Save plot
fig_dir <- "../../tmp/figures/motif_sharing"
output_file <- file.path(fig_dir, "96plex_jaccard_similarity_filtered.pdf")
ggsave(output_file, plot = p, width = 6, height = 4)

cat("Plot saved to:", output_file, "\n")
cat("Number of MGE-host pairs:", nrow(df), "\n")
