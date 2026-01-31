#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(ggplot2)

# This script builds a legend showing MGE types (plasmid, virus) with shapes
# and host phyla as colored circle entries. It saves a PDF legend file.

color_path <- "network_colors.csv"
out_pdf <- "../../tmp/figures/multi_env_linkage/network_99/network_legend.pdf"

# Optional: define explicit colors here. Keys should match category names (plasmid, virus, or phylum strings).
# Example: custom_colors <- c(plasmid = "#A6D854", virus = "#FB9A99", "p__Proteobacteria" = "#66C2A5")
custom_colors <- NULL

if (!file.exists(color_path)) {
  stop(paste("Data file not found:", color_path))
}

df <- read_csv(color_path, show_col_types = FALSE)

# `df` expected columns: Id, Label, type, Polygon, Color
# `type` contains either 'plasmid', 'virus', or taxonomy strings like 'p__Pseudomonadota'

# choose the most frequent color for plasmid and virus (modal color)
plasmid_color <- df %>% filter(type == "plasmid", !is.na(Color)) %>% count(Color, sort = TRUE) %>% slice_head(n = 1) %>% pull(Color)
virus_color <- df %>% filter(type == "virus", !is.na(Color)) %>% count(Color, sort = TRUE) %>% slice_head(n = 1) %>% pull(Color)

# extract phylum entries where `type` starts with 'p__' and pick modal color per phylum
phy_df <- df %>% filter(grepl('^p__', type), !is.na(Color)) %>%
  group_by(type, Color) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(type) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  arrange(type) %>%
  ungroup()
phy_names <- phy_df$type
phy_colors <- phy_df$Color

# Build ordered categories: plasmid, virus, then phyla (if they exist)
categories <- character(0)
if (length(plasmid_color) > 0) categories <- c(categories, "plasmid")
if (length(virus_color) > 0) categories <- c(categories, "virus")
if (length(phy_names) > 0) categories <- c(categories, phy_names)

if (length(categories) == 0) stop("No categories found in color file")

# Shapes: swap plasmid vs hosts: plasmid = circle (16), phyla (hosts) = square (15)
shape_map <- rep(15, length(categories))
names(shape_map) <- categories
# default: phyla/hosts get square (15)
if ("plasmid" %in% categories) shape_map["plasmid"] <- 16
# use triangle for virus
if ("virus" %in% categories) shape_map["virus"] <- 17

# Colors: prefer colors from file; fallback to default palette for any missing
cols <- setNames(rep(NA_character_, length(categories)), categories)
if ("plasmid" %in% categories) cols["plasmid"] <- ifelse(length(plasmid_color)>0, plasmid_color, NA)
if ("virus" %in% categories) cols["virus"] <- ifelse(length(virus_color)>0, virus_color, NA)
for (i in seq_along(phy_names)) {
  cols[phy_names[i]] <- phy_colors[i]
}

library(scales)
missing_cols <- is.na(cols)
if (any(missing_cols)) {
  cols[missing_cols] <- hue_pal()(sum(missing_cols))
}

# overlay custom colors if provided
if (!is.null(custom_colors)) {
  for (nm in names(custom_colors)) if (nm %in% names(cols)) cols[nm] <- custom_colors[[nm]]
}

# Create plotting dataframe
plot_df <- data.frame(category = factor(names(cols), levels = names(cols)), x = 1, y = seq_along(cols))

# Prepare display labels: strip 'p__' from phylum names for clarity
display_labels <- names(cols)
display_labels <- sub('^p__', '', display_labels)

# Build ggplot with both color and shape mapped to the same category name so legends merge
# Make plotted points invisible (alpha=0) so only the legend appears on the right
p <- ggplot(plot_df, aes(x = x, y = y, color = category, shape = category)) +
  geom_point(size = 5, alpha = 0) +
  scale_color_manual(name = "Genome", values = cols, labels = display_labels) +
  scale_shape_manual(name = "Genome", values = shape_map, labels = display_labels) +
  guides(
    color = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 5)),
    shape = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 5))
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
colors_out <- "../../tmp/figures/multi_env_linkage/network_99/network_legend_colors.tsv"
cols_df <- data.frame(category = names(cols), color = as.character(cols), shape = as.integer(shape_map[names(cols)]))
write.table(cols_df, file = colors_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved color mapping to", colors_out, "\n")
