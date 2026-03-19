#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(ggplot2)

# This script builds a legend showing MGE types (plasmid, virus) with shapes
# and host phyla as colored circle entries. It saves a PDF legend file.

color_path <- "network_colors.csv"
out_pdf <- "../../tmp/figures/borg_fig/borg_legend.pdf"

# Optional: define explicit colors here. Keys should match category names (plasmid, virus, or phylum strings).
# Example: custom_colors <- c(plasmid = "#A6D854", virus = "#FB9A99", "p__Proteobacteria" = "#66C2A5")
custom_colors <- NULL

if (!file.exists(color_path)) {
  stop(paste("Data file not found:", color_path))
}

df <- read_csv(color_path, show_col_types = FALSE)

# `df` expected columns: Id, Label, genome, sample, Color
# `genome` contains genome types like 'Mini_Chr', 'Mp', 'Mp_Virus', 'Non-Mp', 'Borg', 'mini_Borg'

# Get unique genome types and their most frequent colors
genome_df <- df %>% 
  filter(!is.na(genome), !is.na(Color)) %>%
  group_by(genome, Color) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(genome) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  arrange(genome) %>%
  ungroup()

genome_names <- genome_df$genome
genome_colors <- genome_df$Color

if (length(genome_names) == 0) stop("No genome types found in color file")

# Build categories from genome types
categories <- genome_names

# All nodes are circles (16)
shape_map <- rep(16, length(categories))
names(shape_map) <- categories

# Assign colors from the data
cols <- setNames(genome_colors, genome_names)

# overlay custom colors if provided
if (!is.null(custom_colors)) {
  for (nm in names(custom_colors)) if (nm %in% names(cols)) cols[nm] <- custom_colors[[nm]]
}

# Create plotting dataframe
plot_df <- data.frame(category = factor(names(cols), levels = names(cols)), x = 1, y = seq_along(cols))

# Use genome names as display labels
display_labels <- names(cols)

# Build ggplot with both color and shape mapped to the same category name so legends merge
# Make plotted points invisible (alpha=0) so only the legend appears on the right
p <- ggplot(plot_df, aes(x = x, y = y, color = category, shape = category)) +
  geom_point(size = 5, alpha = 0) +
  scale_color_manual(name = "Genome", values = cols, labels = display_labels) +
  scale_shape_manual(name = "Genome", values = shape_map, labels = display_labels) +
  guides(
    color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 5)),
    shape = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(alpha = 1, size = 5))
  ) +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# Ensure output directory exists
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# Save legend as PDF
ggsave(out_pdf, plot = p, width = 6, height = length(categories) * 0.25 + 1)
cat("Saved legend to", out_pdf, "\n")

# Also write color mapping to TSV for reproducibility
colors_out <- "../../tmp/figures/borg_fig/borg_legend_colors.tsv"
cols_df <- data.frame(category = names(cols), color = as.character(cols), shape = as.integer(shape_map[names(cols)]))
write.table(cols_df, file = colors_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved color mapping to", colors_out, "\n")
