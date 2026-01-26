#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read data
gc_df <- read_csv("../../tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv", show_col_types = FALSE)

# Rename 'novel' to 'unknown' for consistency and set factor order
gc_df <- gc_df %>%
  mutate(MGE_type_clean = ifelse(MGE_type == "novel", "unknown", MGE_type)) %>%
  mutate(MGE_type_clean = factor(MGE_type_clean, levels = c("plasmid", "virus", "unknown")))

# Count linkages by environment
env_counts <- gc_df %>%
  group_by(environment, MGE_type_clean, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")

# Order environments by total count
env_order <- gc_df %>%
  count(environment, sort = TRUE) %>%
  pull(environment)

env_counts$environment <- factor(env_counts$environment, levels = env_order)

# Count linkages by phylum (top 10)
phylum_counts <- gc_df %>%
  group_by(host_phylum, MGE_type_clean, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")

top_phyla <- gc_df %>%
  count(host_phylum, sort = TRUE) %>%
  head(10) %>%
  pull(host_phylum)

phylum_counts_filtered <- phylum_counts %>%
  filter(host_phylum %in% top_phyla)

phylum_counts_filtered$host_phylum <- factor(phylum_counts_filtered$host_phylum, levels = top_phyla)

# Create color palette
mge_colors <- c("plasmid" = "#0000FF", "virus" = "#FF0000", "unknown" = "#FFA500")

# Plot 1: By environment (without legend)
p1 <- ggplot(env_counts, aes(x = environment, y = count, fill = MGE_type_clean)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = mge_colors, name = "MGE Type") +
  labs(
    x = "Environment",
    y = "Number of Linkages",
    title = "Linkages by Environment"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Plot 2: By phylum (without legend)
p2 <- ggplot(phylum_counts_filtered, aes(x = host_phylum, y = count, fill = MGE_type_clean)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = mge_colors, name = "MGE Type") +
  labs(
    x = "Host Phylum",
    y = "Number of Linkages",
    title = "Linkages by Phylum (Top 10)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Extract legend from a plot with legend
p_legend <- ggplot(env_counts, aes(x = environment, y = count, fill = MGE_type_clean)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = mge_colors, name = "MGE Type") +
  theme_bw() +
  theme(legend.position = "bottom")

# Combine plots with shared legend
library(gridExtra)
library(grid)

# Extract the legend
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p_legend)

# Arrange plots in 2 rows with shared legend at bottom
combined <- grid.arrange(
  arrangeGrob(p1, p2, nrow = 2),
  legend,
  nrow = 2,
  heights = c(10, 1)
)

ggsave("../../tmp/figures/multi_env_linkage/network_99/linkage_counts_combined.pdf", 
       plot = combined, width = 7, height = 8)

# Print summary statistics
cat("\n=== Linkage Count Summary ===\n\n")
cat("By Environment:\n")
print(env_counts %>% pivot_wider(names_from = MGE_type_clean, values_from = count, values_fill = 0))

cat("\nBy Phylum (Top 10):\n")
print(phylum_counts_filtered %>% pivot_wider(names_from = MGE_type_clean, values_from = count, values_fill = 0))

cat("\nPlot saved to ../../tmp/figures/multi_env_linkage/network_99/linkage_counts_combined.pdf\n")
