#!/usr/bin/env Rscript
# Plot detected motif number vs depth for each contig (pure2 96plex subsamples).
# Usage: Rscript plot_motif_vs_depth.R
#   (run from MODIFI/benchmark/base or use full path to script)
# Reads: /home/shuaiw/borg/paper/linkage/pure2/coverage_motif_summary.csv
# Writes: /home/shuaiw/MODIFI/tmp/figures/base_benchmark/

library(ggplot2)
library(tidyr)
library(dplyr)

base_dir <- "/home/shuaiw/MODIFI/tmp/figures/base_benchmark"
data_path <- "/home/shuaiw/borg/paper/linkage/pure2/coverage_motif_summary.csv"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(data_path, stringsAsFactors = FALSE)
# Ensure depth order for plotting (by numeric subsample level)
d <- d %>%
  mutate(
    pct = as.numeric(gsub(".*\\.p([0-9]+)$", "\\1", folder)),
    pct = ifelse(is.na(pct), 0, pct)
  ) %>%
  arrange(contig, pct)

# Remove contigs with very high motif count (outliers that compress the rest)
max_motif_per_contig <- d %>% group_by(contig) %>% summarise(max_motif = max(motif_count, na.rm = TRUE), .groups = "drop")
contigs_exclude <- max_motif_per_contig %>% filter(max_motif > 120) %>% pull(contig)
if (length(contigs_exclude) > 0) {
  message("Excluding high-motif contigs: ", paste(contigs_exclude, collapse = ", "))
  d <- d %>% filter(!contig %in% contigs_exclude)
}

# Restrict to depth 0-100 for plotting
d_plot <- d %>% filter(depth >= 0, depth <= 100)

# Base = max motif count at depth > 50 per contig; remove contigs with no >50x or base 0
base_at_high <- d_plot %>% filter(depth > 50) %>% group_by(contig) %>% summarise(base = max(motif_count, na.rm = TRUE), .groups = "drop")
contigs_with_high <- base_at_high %>% filter(base > 0) %>% pull(contig)
d_plot <- d_plot %>% filter(contig %in% contigs_with_high) %>% left_join(base_at_high, by = "contig")
d_plot <- d_plot %>% mutate(motif_pct = 100 * motif_count / base)
if (length(contigs_with_high) == 0) stop("No contigs with depth > 50 and motif count > 0")
message("Contigs with >50x baseline: ", length(contigs_with_high), " (removed ", n_distinct(d$contig) - length(contigs_with_high), " without >50x)")

# Order contigs by base (descending) for legend/ordering
contig_order <- base_at_high %>% filter(contig %in% contigs_with_high) %>% arrange(desc(base)) %>% pull(contig)
d_plot <- d_plot %>% mutate(contig = factor(contig, levels = rev(contig_order)))

# ---- 1) Line plot: motif count as % of >50x baseline (depth 0-100) ----
p1 <- ggplot(d_plot, aes(x = depth, y = motif_pct, color = contig)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  geom_point(size = 1.4, alpha = 0.9) +
  scale_x_continuous("Depth", limits = c(0, 100), breaks = seq(0, 100, by = 10), expand = c(0.02, 0)) +
  scale_y_continuous("Motif count (% of >50x)", limits = c(0, NA), expand = c(0, 2)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.15, colour = "grey92"),
    plot.title = element_text(size = 13, face = "bold")
  ) +
  ggtitle("Motif count vs depth (% of >50x baseline)")

out1 <- file.path(base_dir, "motif_count_vs_depth_all.pdf")
ggsave(out1, p1, width = 8, height = 3.5)
message("Saved: ", out1)

out1_png <- file.path(base_dir, "motif_count_vs_depth_all.png")
ggsave(out1_png, p1, width = 8, height = 3.5, dpi = 200)
message("Saved: ", out1_png)

# ---- 1b) Heatmap: contigs x depth, fill = motif % of >50x ----
p_heat <- ggplot(d_plot, aes(x = depth, y = contig, fill = motif_pct)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_viridis_c("Motif (% of >50x)", option = "plasma", begin = 0.1, end = 0.95, limits = c(0, 100)) +
  scale_x_continuous("Depth", limits = c(0, 100), breaks = seq(0, 100, by = 5), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold")
  ) +
  ggtitle("Motif count vs depth (% of >50x baseline)")

out_heat <- file.path(base_dir, "motif_count_vs_depth_heatmap.pdf")
ggsave(out_heat, p_heat, width = 6.5, height = max(4, nlevels(d_plot$contig) * 0.2))
message("Saved: ", out_heat)

out_heat_png <- file.path(base_dir, "motif_count_vs_depth_heatmap.png")
ggsave(out_heat_png, p_heat, width = 6.5, height = max(4, nlevels(d_plot$contig) * 0.2), dpi = 200)
message("Saved: ", out_heat_png)

# ---- 2) Faceted: one panel per contig (% of >50x, depth 0-100) ----
p2 <- ggplot(d_plot, aes(x = depth, y = motif_pct)) +
  geom_line(color = "steelblue", linewidth = 0.7) +
  geom_point(color = "steelblue", size = 1.5) +
  facet_wrap(~ contig, scales = "free_y", ncol = 4) +
  scale_x_continuous("Depth", limits = c(0, 60), breaks = seq(0, 60, by = 5)) +
  scale_y_continuous("Motif count (% of >50x)") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 9),
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.15, colour = "grey92"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ggtitle("Motif count vs depth by contig (% of >50x)")

out2 <- file.path(base_dir, "motif_count_vs_depth_by_contig.pdf")
ggsave(out2, p2, width = 8, height = 6)
message("Saved: ", out2)

out2_png <- file.path(base_dir, "motif_count_vs_depth_by_contig.png")
ggsave(out2_png, p2, width = 8, height = 6, dpi = 150)
message("Saved: ", out2_png)

# ---- 3) Summary: mean motif % by depth bin (0-10, 10-20, ...) with error bars (std) ----
depth_bin_width <- 10
d_summary <- d_plot %>%
  mutate(depth_bin = depth_bin_width * floor(depth / depth_bin_width)) %>%
  group_by(depth_bin) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    mean_pct = mean(motif_pct, na.rm = TRUE),
    sd_pct = sd(motif_pct, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    sd_pct = ifelse(is.na(sd_pct) | n < 2, 0, sd_pct),
    ymin = pmax(0, mean_pct - sd_pct),
    ymax = mean_pct + sd_pct
  )

# Save mean plot data to figure folder
out_data <- file.path(base_dir, "motif_count_vs_depth_mean_data.csv")
write.csv(d_summary, out_data, row.names = FALSE)
message("Saved: ", out_data)

p3 <- ggplot(d_summary, aes(x = mean_depth, y = mean_pct)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), color = "darkred", linewidth = 0.6, width = 2) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_point(color = "darkred", size = 3) +
  scale_x_continuous("Depth (x)", limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous("Motif detection recovery rate (%)", limits = c(0, 110), breaks = seq(0, 100, by = 20)) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
    panel.grid.minor = element_line(linewidth = 0.15, colour = "grey92")
  ) 

out3 <- file.path(base_dir, "motif_count_vs_depth_mean.pdf")
ggsave(out3, p3, width = 5, height = 3.5)
message("Saved: ", out3)

message("Done. Figures in: ", base_dir)
