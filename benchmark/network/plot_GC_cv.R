#!/usr/bin/env Rscript
# Converted from plot_linkage_data.py plot_gc() (MGE vs host GC, coverage, cosine similarity)

library(ggplot2)
library(patchwork)
library(dplyr)

# Environment color map (same as plot_motif_MTase_corr.R)
env_colors <- c(
  mock = "#e41a1c",
  mice_gut = "#377eb8",
  sugarcane = "#4daf4a",
  infant_gut = "#984ea3",
  cow_rumen = "#ff7f00",
  cow_bioreactor = "#ffff33",
  adult_gut = "#a65628",
  ocean = "#f781bf",
  soil = "#999999"
)

# Path to CSV produced by plot_linkage_data.py and output directory
# Default: same directory as script; override as needed
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  paper_fig_dir <- args[1]
} else {
  paper_fig_dir <- "../../tmp/figures/multi_env_linkage/"
}
csv_path <- paste0(paper_fig_dir, "/network_99/mge_host_gc_cov.csv")

if (!file.exists(csv_path)) {
  stop("CSV not found: ", csv_path, "\nGenerate it by running plot_linkage_data.py first.")
}

df <- read.csv(csv_path, stringsAsFactors = FALSE)
required <- c("sample", "MGE_gc", "host_gc", "MGE_cov", "host_cov", "cos_sim", "environment", "MGE_type")
missing <- setdiff(required, names(df))
if (length(missing) > 0) {
  stop("CSV missing columns: ", paste(missing, collapse = ", "))
}
df <- df %>% arrange(sample)

# MGE_type shape mapping (plasmid=circle, virus=X, novel=triangle); use only levels present
shape_map_all <- c(plasmid = 16, virus = 4, novel = 17)
mge_levels <- intersect(names(shape_map_all), unique(as.character(df$MGE_type[!is.na(df$MGE_type)])))
if (length(mge_levels) == 0) mge_levels <- unique(as.character(df$MGE_type))
if (length(mge_levels) == 0) mge_levels <- "plasmid"
df$MGE_type <- factor(df$MGE_type, levels = mge_levels)
shape_map <- shape_map_all[names(shape_map_all) %in% mge_levels]

# Extend env colors for any environments not in the map
env_levels <- unique(as.character(df$environment))
extra_envs <- setdiff(env_levels, names(env_colors))
if (length(extra_envs) > 0) {
  env_colors <- c(env_colors, setNames(rep("#999999", length(extra_envs)), extra_envs))
}

# ---- Panel 1: MGE_gc vs host_gc ----
valid_gc <- df %>% filter(!is.na(MGE_gc), !is.na(host_gc))
if (nrow(valid_gc) > 1) {
  ct_gc <- cor.test(valid_gc$MGE_gc, valid_gc$host_gc, method = "pearson", exact = FALSE)
  r_gc <- ct_gc$estimate
  p_gc <- ct_gc$p.value
} else {
  r_gc <- NA
  p_gc <- NA
}

gc_vals <- c(df$MGE_gc, df$host_gc)
gc_vals <- gc_vals[!is.na(gc_vals)]
if (length(gc_vals) == 0) {
  min_gc <- 0
  max_gc <- 1
} else {
  min_gc <- min(gc_vals)
  max_gc <- max(gc_vals)
}

p1 <- ggplot(df, aes(x = MGE_gc, y = host_gc, color = environment, shape = MGE_type)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8, alpha = 0.15, inherit.aes = FALSE, aes(x = MGE_gc, y = host_gc), data = df) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = env_colors, name = "Habitat") +
  scale_shape_manual(values = shape_map, na.value = 16, drop = TRUE) +
  scale_x_continuous(limits = c(min_gc, max_gc)) +
  scale_y_continuous(limits = c(min_gc, max_gc)) +
  labs(x = "MGE GC content", y = "Host GC content") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(5, 5, 5, 5)
  )
if (!is.na(r_gc)) {
  p1 <- p1 + annotate("text", x = min_gc, y = max_gc,
                      label = sprintf("r = %.2f\np = %.2g", r_gc, p_gc),
                      hjust = 0, vjust = 1, size = 3)
}

# ---- Panel 2: MGE_cov vs host_cov (log scale) ----
valid_cov <- df %>% filter(MGE_cov > 0, host_cov > 0)
if (nrow(valid_cov) > 1) {
  ct_cov <- cor.test(log10(valid_cov$MGE_cov), log10(valid_cov$host_cov), method = "pearson", exact = FALSE)
  r_cov <- ct_cov$estimate
  p_cov <- ct_cov$p.value
} else {
  r_cov <- NA
  p_cov <- NA
}

df_cov <- df %>% filter(MGE_cov > 0, host_cov > 0)
# Regression line and CI band in log-log space (straight on log-scale plot)
if (nrow(df_cov) > 1) {
  fit_cov <- lm(log10(host_cov) ~ log10(MGE_cov), data = df_cov)
  xseq <- seq(min(df_cov$MGE_cov), max(df_cov$MGE_cov), length.out = 100)
  pred_cov <- predict(fit_cov, newdata = data.frame(MGE_cov = xseq), interval = "confidence")
  line_cov <- data.frame(
    MGE_cov = xseq,
    host_cov = 10^pred_cov[, "fit"],
    host_cov_lwr = 10^pred_cov[, "lwr"],
    host_cov_upr = 10^pred_cov[, "upr"]
  )
} else {
  line_cov <- data.frame(MGE_cov = numeric(0), host_cov = numeric(0), host_cov_lwr = numeric(0), host_cov_upr = numeric(0))
}
p2 <- ggplot(df_cov, aes(x = MGE_cov, y = host_cov, color = environment, shape = MGE_type)) +
  geom_ribbon(data = line_cov, aes(x = MGE_cov, ymin = host_cov_lwr, ymax = host_cov_upr), fill = "black", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = line_cov, aes(x = MGE_cov, y = host_cov), color = "black", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = env_colors, name = "Habitat") +
  scale_shape_manual(values = shape_map, na.value = 16) +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "MGE Coverage (log scale)", y = "Host Coverage (log scale)") +
  theme_minimal() +
  theme(legend.position = "none", plot.margin = margin(5, 5, 5, 5))
if (!is.na(r_cov) && nrow(valid_cov) > 0) {
  p2 <- p2 + annotate("text", x = min(valid_cov$MGE_cov), y = max(valid_cov$host_cov),
                      label = sprintf("r = %.2f\np = %.2g", r_cov, p_cov),
                      hjust = 0, vjust = 1, size = 3)
}

# ---- Panel 3: cos_sim by environment (boxplot) ----
env_medians <- df %>%
  group_by(environment) %>%
  summarise(med = median(cos_sim, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med))
env_order <- env_medians$environment
if (length(env_order) == 0) env_order <- unique(df$environment)
df$environment_ordered <- factor(df$environment, levels = env_order)

p3 <- ggplot(df, aes(x = environment_ordered, y = cos_sim, fill = environment)) +
  geom_boxplot(outlier.size = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = env_colors, name = "Habitat") +
  labs(x = "Habitat", y = "Cosine Similarity") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.margin = margin(5, 5, 5, 5)
  )

# ---- Panel 4: empty ----
p4 <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
  geom_blank() +
  theme_void() +
  theme(panel.border = element_blank())

# ---- Combine 2x2 and save ----
combined <- (p1 | p2) / (p3 | p4) +
  plot_layout(heights = c(1, 1), widths = c(1, 1))

out_path <- file.path("../../tmp/figures/multi_env_linkage/network_99/mge_host_gc_content.pdf")
ggsave(out_path, combined, width = 14, height = 10, dpi = 300)
cat("Saved:", out_path, "\n")
