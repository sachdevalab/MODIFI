#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# Read data
data_file <- "/home/shuaiw/mGlu/tmp/figures/motif_sharing/jaccard_similarity_96plex.csv"
df <- read_csv(data_file, show_col_types = FALSE)

# Reshape data to long format for grouped bars
df_long <- df %>%
  select(MGE, host, jaccard_similarity, jaccard_similarity_filtered, 
         total_motif_num, total_motif_num_filtered) %>%
  pivot_longer(cols = c(jaccard_similarity, jaccard_similarity_filtered),
               names_to = "similarity_type",
               values_to = "similarity_value") %>%
  mutate(
    motif_count = ifelse(similarity_type == "jaccard_similarity", 
                         total_motif_num, total_motif_num_filtered),
    type_label = factor(ifelse(similarity_type == "jaccard_similarity",
                       "Original", "Filtered"),
                       levels = c("Original", "Filtered"))
  )

# Reorder MGE to put S_enterica_LT2_2 first
mge_order <- unique(df_long$MGE)
mge_order <- c("S_enterica_LT2_2", setdiff(mge_order, "S_enterica_LT2_2"))
df_long$MGE <- factor(df_long$MGE, levels = mge_order)

# Create the plot
p <- ggplot(df_long, aes(x = MGE, y = similarity_value, fill = type_label)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = motif_count), 
            position = position_dodge(width = 0.8),
            vjust = 1.5, size = 3, color = "white") +
  scale_fill_manual(values = c("Original" = "#8B0000", "Filtered" = "#5F9EA0"),
                    name = "Motif filtering") +
  labs(x = "",
       y = "Jaccard Similarity",
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.title = element_text(size = 12, ),
        plot.title = element_text(size = 14, hjust = 0.5),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)) +
  ylim(0, 1)

# Save plot
fig_dir <- "../../tmp/figures/motif_sharing"
output_file <- file.path(fig_dir, "96plex_jaccard_comparison.pdf")
ggsave(output_file, plot = p, width = 8, height = 4)

cat("Plot saved to:", output_file, "\n")
cat("Number of MGE-host pairs:", nrow(df), "\n")
