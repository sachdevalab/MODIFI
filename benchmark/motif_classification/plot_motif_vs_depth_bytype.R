#!/usr/bin/env Rscript
# Panel d (PB24) split by modification type: 6mA / 4mC / 5mC / unclassified.
# Fork of plot_motif_vs_depth.R panel p3 (recovery rate vs depth). Original untouched.
#
# Recovery per type = motif_count normalised to that (contig, type)'s own >50x baseline,
# then averaged across contigs per depth bin (mean +/- sd), exactly as the original panel d
# but stratified by modification type. 5mC has 0 detected motifs at any depth, so it is drawn
# as a flat 0% line (labelled "5mC (0 detected)") to make the PacBio 5mC blind spot explicit.
#
# Reads : /home/shuaiw/borg/revision/motif_class/coverage_motif_summary_bytype.csv
# Writes: /home/shuaiw/MODIFI/tmp/rev_figs/motif_classification/
#           panel_d_by_modtype.pdf/.png  +  panel_d_by_modtype_sourcedata.csv

library(ggplot2)
library(dplyr)

data_path <- "/home/shuaiw/borg/revision/motif_class/coverage_motif_summary_bytype.csv"
out_dir   <- "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

types <- c("6mA", "4mC", "5mC", "unclassified")
type_cols <- c("6mA" = "#C0392B", "4mC" = "#2E86C1",
               "5mC" = "#7F7F7F", "unclassified" = "#E59400")

d <- read.csv(data_path, stringsAsFactors = FALSE)
d <- d %>% mutate(pct = as.numeric(gsub(".*\\.p([0-9]+)$", "\\1", folder)),
                  pct = ifelse(is.na(pct), 0, pct))

# Match the original panel d contig set: drop contigs whose TOTAL motif count (summed over
# types, max over folders) exceeds 120 (compressing outliers).
tot_per_fc <- d %>% group_by(folder, contig) %>%
  summarise(tot = sum(motif_count), .groups = "drop")
excl <- tot_per_fc %>% group_by(contig) %>% summarise(mx = max(tot), .groups = "drop") %>%
  filter(mx > 120) %>% pull(contig)
if (length(excl) > 0) message("Excluding high-motif contigs: ", paste(excl, collapse = ", "))
d <- d %>% filter(!contig %in% excl)

# Depth window 0-100
d_plot <- d %>% filter(depth >= 0, depth <= 100)

# Per (contig, type) baseline = max count at depth > 50; keep only those with base > 0.
base_ct <- d_plot %>% filter(depth > 50) %>%
  group_by(contig, mod_type) %>%
  summarise(base = max(motif_count, na.rm = TRUE), .groups = "drop") %>%
  filter(base > 0)

d_rec <- d_plot %>% inner_join(base_ct, by = c("contig", "mod_type")) %>%
  mutate(motif_pct = 100 * motif_count / base)

# Report how many contigs contribute to each type
ct_n <- d_rec %>% distinct(contig, mod_type) %>% count(mod_type)
message("Contigs contributing per type:")
for (i in seq_len(nrow(ct_n))) message("  ", ct_n$mod_type[i], ": ", ct_n$n[i])

# Aggregate: mean +/- sd of recovery across contigs, per 10x depth bin and type
bw <- 10
d_sum <- d_rec %>%
  mutate(depth_bin = bw * floor(depth / bw)) %>%
  group_by(mod_type, depth_bin) %>%
  summarise(mean_depth = mean(depth, na.rm = TRUE),
            mean_pct = mean(motif_pct, na.rm = TRUE),
            sd_pct = sd(motif_pct, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(sd_pct = ifelse(is.na(sd_pct) | n < 2, 0, sd_pct),
         ymin = pmax(0, mean_pct - sd_pct),
         ymax = pmin(110, mean_pct + sd_pct))

# 5mC has 0 detected motifs at any depth, so no 5mC curve is drawn (dropped per request).
plot_types <- c("6mA", "4mC")

# Legend labels carry the contributing-contig count
lab_n <- setNames(rep("", length(plot_types)), plot_types)
for (t in plot_types) {
  k <- ct_n$n[ct_n$mod_type == t]
  lab_n[t] <- if (length(k) == 0) t else paste0(t, " (n=", k, ")")
}
d_sum <- d_sum %>% filter(mod_type %in% plot_types) %>%
  mutate(mod_type = factor(mod_type, levels = plot_types))

# Fixed horizontal offset per type so the two series' points/error bars don't overlap.
x_off <- c("6mA" = -1.3, "4mC" = 1.3)
d_sum <- d_sum %>% mutate(x_plot = mean_depth + x_off[as.character(mod_type)])

write.csv(d_sum, file.path(out_dir, "panel_d_by_modtype_sourcedata.csv"), row.names = FALSE)

p <- ggplot(d_sum, aes(x = x_plot, y = mean_pct, color = mod_type)) +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), linewidth = 0.5, width = 2, alpha = 0.9) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  scale_color_manual("Modification type", values = type_cols,
                     breaks = plot_types, labels = lab_n[plot_types]) +
  scale_x_continuous("Depth (x)", limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  scale_y_continuous("Motif detection recovery rate (%)",
                     limits = c(0, 110), breaks = seq(0, 100, by = 20)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
        panel.grid.minor = element_line(linewidth = 0.15, colour = "grey92"),
        legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))

ggsave(file.path(out_dir, "panel_d_by_modtype.pdf"), p, width = 6, height = 3.8)
ggsave(file.path(out_dir, "panel_d_by_modtype.png"), p, width = 6, height = 3.8, dpi = 200)
message("Saved: ", file.path(out_dir, "panel_d_by_modtype.pdf"))
