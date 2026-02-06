#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(patchwork)

# Read the melted matrix
melt_df <- read_csv("../../tmp/figures/multi_env_linkage/network_99/largest_component_binary_matrix_melted.csv", show_col_types = FALSE)

# Convert melted data back to binary matrix
binary_matrix <- melt_df %>%
  select(MGE_node, host_node, present) %>%
  pivot_wider(names_from = host_node, values_from = present, values_fill = 0) %>%
  column_to_rownames("MGE_node")

# Create MGE type and label mappings
mge_info <- melt_df %>%
  select(MGE_node, MGE_type, MGE_label) %>%
  distinct()

host_info <- melt_df %>%
  select(host_node, host_label, host_species) %>%
  distinct()

# Separate virus and non-virus rows
virus_rows <- mge_info %>% filter(MGE_type == "virus") %>% pull(MGE_node)
non_virus_rows <- mge_info %>% filter(MGE_type != "virus") %>% pull(MGE_node)

# Filter to rows that exist in the binary matrix
virus_rows <- intersect(virus_rows, rownames(binary_matrix))
non_virus_rows <- intersect(non_virus_rows, rownames(binary_matrix))

# Perform hierarchical clustering on non-virus rows
if (length(non_virus_rows) > 1) {
  sub_nv <- binary_matrix[non_virus_rows, , drop = FALSE]
  
  # Cluster rows (MGEs) using Jaccard distance
  row_dist <- dist(sub_nv, method = "binary")  # binary method is equivalent to Jaccard
  row_clust <- hclust(row_dist, method = "average")
  row_order <- row_clust$order
} else {
  row_order <- seq_along(non_virus_rows)
}

# Cluster columns (hosts) based on non-virus rows
if (length(non_virus_rows) > 0 && ncol(binary_matrix) > 1) {
  col_dist <- dist(t(binary_matrix[non_virus_rows, , drop = FALSE]), method = "binary")
  col_clust <- hclust(col_dist, method = "average")
  col_order <- col_clust$order
} else {
  col_order <- seq_len(ncol(binary_matrix))
}

# Reorder matrix
ordered_nonvirus <- binary_matrix[non_virus_rows[row_order], colnames(binary_matrix)[col_order], drop = FALSE]

if (length(virus_rows) > 0) {
  ordered_virus <- binary_matrix[virus_rows, colnames(binary_matrix)[col_order], drop = FALSE]
  binary_ordered <- rbind(ordered_nonvirus, ordered_virus)
} else {
  binary_ordered <- ordered_nonvirus
}

# Create host labels with species information
host_labels_df <- host_info %>%
  filter(host_node %in% colnames(binary_ordered)) %>%
  mutate(label = paste0(host_species, " (", host_label, ")"))

host_labels <- setNames(host_labels_df$label, host_labels_df$host_node)
ordered_col_labels <- host_labels[colnames(binary_ordered)]

# Prepare data for ggplot
plot_data <- binary_ordered %>%
  as.data.frame() %>%
  rownames_to_column("MGE_node") %>%
  pivot_longer(-MGE_node, names_to = "host_node", values_to = "present") %>%
  left_join(mge_info, by = "MGE_node") %>%
  mutate(
    MGE_node = factor(MGE_node, levels = rownames(binary_ordered)),
    host_node = factor(host_node, levels = colnames(binary_ordered))
  )

# Create color mapping for MGE types
mge_colors <- setNames(
  c("blue", "red", "orange", "black"),
  c("plasmid", "virus", "novel", "other")
)

# Get row label colors
row_label_colors <- mge_info %>%
  filter(MGE_node %in% rownames(binary_ordered)) %>%
  mutate(
    MGE_node = factor(MGE_node, levels = rownames(binary_ordered)),
    color = case_when(
      MGE_type == "plasmid" ~ "#F8766D",
      MGE_type == "virus" ~ "#00BFC4",
      MGE_type == "novel" ~ "orange",
      TRUE ~ "black"
    )
  ) %>%
  arrange(MGE_node)


# Calculate column sums (number of links per host)
col_sums <- colSums(binary_ordered)
col_sums_df <- data.frame(
  host_node = factor(names(col_sums), levels = colnames(binary_ordered)),
  count = col_sums
)

# Calculate row sums (number of links per MGE)
row_sums <- rowSums(binary_ordered)
row_sums_df <- data.frame(
  MGE_node = factor(names(row_sums), levels = rownames(binary_ordered)),
  count = row_sums
)

# Prepare MGE type annotation data
mge_type_df <- mge_info %>%
  filter(MGE_node %in% rownames(binary_ordered)) %>%
  mutate(MGE_node = factor(MGE_node, levels = rownames(binary_ordered))) %>%
  arrange(MGE_node)

# Define colors for MGE types (matching the provided legend)
type_colors <- c("plasmid" = "#F8766D", "virus" = "#00BFC4", "novel" = "orange")

# Create x-axis labels plot (MGE labels, rotated)
p_xlabel <- ggplot(plot_data %>% distinct(MGE_node, .keep_all = TRUE), 
                   aes(x = MGE_node, y = 1)) +
  geom_text(aes(label = MGE_label, color = MGE_type), angle = 90, hjust = 1, vjust = 0.5, size = 2.2) +
  scale_color_manual(values = c("plasmid" = "blue", "virus" = "red", "novel" = "orange")) +
  scale_x_discrete(drop = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 5, 5)
  )

# Create heatmap (rotated: hosts on y-axis, MGEs on x-axis)
p <- ggplot(plot_data, aes(x = MGE_node, y = host_node, fill = factor(present))) +
  geom_tile(color = "gray90", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "white", "1" = "black"), guide = "none") +
  scale_y_discrete(labels = ordered_col_labels) +
  scale_x_discrete(labels = plot_data %>% distinct(MGE_node, .keep_all = TRUE) %>% pull(MGE_label)) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6.5, 
                                color = row_label_colors$color),
    axis.text.y = element_text(size = 7),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    plot.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 5, 8)
  )

# Create top bar plot (MGE counts)
p_top <- ggplot(row_sums_df, aes(x = MGE_node, y = count)) +
  geom_bar(stat = "identity", fill = "gray40", width = 1) +
  labs(y = "") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 2, 5)
  )

# Create left bar plot (host counts)
p_left <- ggplot(col_sums_df, aes(x = count, y = host_node)) +
  geom_bar(stat = "identity", fill = "gray40", width = 1) +
  labs(x = "") +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 5, 2, 2)
  )

# Create annotation bar (MGE type, horizontal)
p_annot <- ggplot(mge_type_df, aes(x = MGE_node, y = 1, fill = MGE_type)) +
  geom_tile(color = NA) +
  scale_fill_manual(values = type_colors, name = "MGE Type") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 5, 2, 5)
  )

# Combine plots using patchwork
combined_plot <- p_top + plot_spacer() + 
                 p + p_left +
                 plot_layout(
                   ncol = 2, nrow = 2,
                   widths = c(5, 0.5),
                   heights = c(0.8, 3)
                 )

# Save plot
ggsave("../../tmp/figures/multi_env_linkage/network_99/largest_component_binary_heatmap_from_R.pdf",
       plot = combined_plot, width = 9, height = 6, dpi = 300)

cat("Saved clustered heatmap to largest_component_binary_heatmap_from_R.pdf\n")
