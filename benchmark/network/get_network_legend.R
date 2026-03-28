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

# MGE categories (plasmid, virus) and phylum categories separately
mge_categories <- character(0)
if (length(plasmid_color) > 0) mge_categories <- c(mge_categories, "plasmid")
if (length(virus_color) > 0) mge_categories <- c(mge_categories, "virus")
phy_categories <- phy_names

if (length(mge_categories) == 0 && length(phy_categories) == 0) stop("No categories found in color file")

# MGE: shapes and colors
mge_shape_map <- setNames(c(16, 17)[seq_along(mge_categories)], mge_categories)  # circle, triangle
mge_cols <- setNames(character(length(mge_categories)), mge_categories)
if ("plasmid" %in% mge_categories) mge_cols["plasmid"] <- plasmid_color
if ("virus" %in% mge_categories) mge_cols["virus"] <- virus_color

# Phyla: colors (all squares in second legend)
phy_cols <- setNames(phy_colors, phy_categories)

library(scales)
if (any(is.na(mge_cols))) mge_cols[is.na(mge_cols)] <- hue_pal()(sum(is.na(mge_cols)))
if (!is.null(custom_colors)) {
  for (nm in names(custom_colors)) if (nm %in% names(mge_cols)) mge_cols[nm] <- custom_colors[[nm]]
  for (nm in names(custom_colors)) if (nm %in% names(phy_cols)) phy_cols[nm] <- custom_colors[[nm]]
}

# Data for legend only: one row per MGE, then one row per phylum (invisible points)
mge_df <- data.frame(category = factor(mge_categories, levels = mge_categories), x = 1, y = seq_along(mge_categories))
phy_df <- data.frame(category = factor(phy_categories, levels = phy_categories), x = 1, y = seq_along(phy_categories))

mge_labels <- sub('^p__', '', mge_categories)
phy_labels <- sub('^p__', '', phy_categories)

# Two legends: MGE (color + shape), Host phylum (fill)
p <- ggplot() +
  geom_point(data = mge_df, aes(x = x, y = y, color = category, shape = category), size = 5, alpha = 0) +
  scale_color_manual(name = "MGE", values = mge_cols, labels = mge_labels,
                     guide = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 5))) +
  scale_shape_manual(name = "MGE", values = mge_shape_map, labels = mge_labels,
                     guide = guide_legend(ncol = 2, override.aes = list(alpha = 1, size = 5)))

if (length(phy_categories) > 0) {
  p <- p +
    geom_point(data = phy_df, aes(x = x, y = y, fill = category), shape = 22, size = 5, alpha = 0) +
    scale_fill_manual(name = "Host phylum", values = phy_cols, labels = phy_labels,
                      guide = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(alpha = 1, size = 5, shape = 22)))
}

p <- p +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, family = "Arial"),
        legend.text = element_text(size = 10, family = "Arial"))

# Ensure output directory exists
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

# Save legend as PDF
ggsave(out_pdf, plot = p, width = 6, height = (length(mge_categories) + length(phy_categories)) * 0.25 + 1, device = cairo_pdf)
cat("Saved legend to", out_pdf, "\n")

# Also write color mapping to TSV for reproducibility
colors_out <- "../../tmp/figures/multi_env_linkage/network_99/network_legend_colors.tsv"
mge_df_out <- if (length(mge_categories) > 0)
  data.frame(legend = "MGE", category = names(mge_cols), color = as.character(mge_cols), shape = as.integer(mge_shape_map[names(mge_cols)])) else data.frame(legend = character(), category = character(), color = character(), shape = integer())
phy_df_out <- if (length(phy_categories) > 0)
  data.frame(legend = "Host phylum", category = names(phy_cols), color = as.character(phy_cols), shape = 22L) else data.frame(legend = character(), category = character(), color = character(), shape = integer())
cols_df <- rbind(mge_df_out, phy_df_out)
write.table(cols_df, file = colors_out, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Saved color mapping to", colors_out, "\n")
