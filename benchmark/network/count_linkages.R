#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read data
gc_df <- read_csv("../../tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv", show_col_types = FALSE)

gc_df <- gc_df %>%
  mutate(MGE_type = factor(MGE_type, levels = c("plasmid", "virus")))


# Count linkages per sample by environment and MGE type
env_sample_counts <- gc_df %>%
  group_by(environment, sample, MGE_type, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")

# Only retain environments with >=5 samples (for both MGE types combined)
env_sample_counts <- env_sample_counts %>%
  group_by(environment) %>%
  filter(n_distinct(sample) >= 3) %>%
  ungroup()

# Order environments by total count (for consistent x-axis order)
env_order <- gc_df %>%
  count(environment, sort = TRUE) %>%
  pull(environment)
env_sample_counts$environment <- factor(env_sample_counts$environment, levels = env_order)

# Count linkages per sample by phylum and MGE type (top 10 phyla)
top_phyla <- gc_df %>%
  count(host_phylum, sort = TRUE) %>%
  head(10) %>%
  pull(host_phylum)

phylum_sample_counts <- gc_df %>%
  filter(host_phylum %in% top_phyla) %>%
  mutate(host_phylum = gsub("^p__", "", host_phylum)) %>%
  group_by(host_phylum, sample, MGE_type, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")
phylum_sample_counts$host_phylum <- factor(phylum_sample_counts$host_phylum, levels = gsub("^p__", "", top_phyla))

# Only retain phyla with >=5 samples (for both MGE types combined)
phylum_sample_counts <- phylum_sample_counts %>%
  group_by(host_phylum) %>%
  filter(n_distinct(sample) >= 5) %>%
  ungroup()


# Load ggpubr for significance annotation
if (!requireNamespace("ggpubr", quietly = TRUE)) {
  install.packages("ggpubr", repos = "https://cloud.r-project.org")
}
library(ggpubr)
# Ensure symnum is available (from stats)
if (!exists("symnum")) symnum <- get("symnum", envir = asNamespace("stats"))


# Boxplot 1: Linkages per sample by environment (filtered, with significance)
# Precompute p-values for each environment
env_pvals <- env_sample_counts %>%
  group_by(environment) %>%
  filter(length(unique(MGE_type)) == 2) %>%
  summarise(p = tryCatch(t.test(count ~ MGE_type)$p.value, error = function(e) NA_real_)) %>%
  mutate(star = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "")))

# Get y position for annotation
env_ypos <- env_sample_counts %>% group_by(environment) %>% summarise(y = max(count, na.rm = TRUE) * 1.1)
env_pvals <- left_join(env_pvals, env_ypos, by = "environment")

p1 <- ggplot(env_sample_counts, aes(x = environment, y = count, fill = MGE_type)) +
  geom_boxplot(aes(fill = MGE_type), outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 1, color = "black", size = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 1.2, alpha = 1, show.legend = FALSE) +
  scale_fill_discrete(name = "MGE Type") +
  scale_color_discrete(guide = "none") +
  labs(
    x = "Environment",
    y = "Linkages per Sample",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  geom_text(data = env_pvals, aes(x = environment, y = y, label = star), inherit.aes = FALSE, vjust = 0, size = 6)

ggsave("../../tmp/figures/multi_env_linkage/network_99/linkage_counts_env_boxplot.pdf", 
       plot = p1, width = 3.5, height = 6)

# Boxplot 2: Linkages per sample by phylum (filtered, with significance)
# Precompute p-values for each phylum
phylum_pvals <- phylum_sample_counts %>%
  group_by(host_phylum) %>%
  filter(length(unique(MGE_type)) == 2) %>%
  summarise(p = tryCatch(t.test(count ~ MGE_type)$p.value, error = function(e) NA_real_)) %>%
  mutate(star = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "")))

# Get y position for annotation
phylum_ypos <- phylum_sample_counts %>% group_by(host_phylum) %>% summarise(y = max(count, na.rm = TRUE) * 1.1)
phylum_pvals <- left_join(phylum_pvals, phylum_ypos, by = "host_phylum")

p2 <- ggplot(phylum_sample_counts, aes(x = host_phylum, y = count, fill = MGE_type)) +
  geom_boxplot(aes(fill = MGE_type), outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 1, color = "black", size = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 1.2, alpha = 1, show.legend = FALSE) +
  scale_fill_discrete(name = "MGE Type") +
  scale_color_discrete(guide = "none") +
  labs(
    x = "Host Phylum",
    y = "Linkages per Sample",
    title = ""
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  geom_text(data = phylum_pvals, aes(x = host_phylum, y = y, label = star), inherit.aes = FALSE, vjust = 0, size = 6)

ggsave("../../tmp/figures/multi_env_linkage/network_99/linkage_counts_phylum_boxplot.pdf", 
       plot = p2, width = 7, height = 5)


# T-tests for environment
cat("\n=== T-tests: Plasmid vs Virus by Environment (≥5 samples) ===\n")
for (env in unique(env_sample_counts$environment)) {
  dat <- env_sample_counts %>% filter(environment == env)
  if (length(unique(dat$MGE_type)) == 2) {
    t_res <- try(t.test(count ~ MGE_type, data = dat), silent = TRUE)
    if (!inherits(t_res, "try-error")) {
      cat(sprintf("%s: p = %.4g\n", as.character(env), t_res$p.value))
    }
  }
}

# T-tests for phylum
cat("\n=== T-tests: Plasmid vs Virus by Phylum (≥5 samples) ===\n")
for (phy in unique(phylum_sample_counts$host_phylum)) {
  dat <- phylum_sample_counts %>% filter(host_phylum == phy)
  if (length(unique(dat$MGE_type)) == 2) {
    t_res <- try(t.test(count ~ MGE_type, data = dat), silent = TRUE)
    if (!inherits(t_res, "try-error")) {
      cat(sprintf("%s: p = %.4g\n", as.character(phy), t_res$p.value))
    }
  }
}

# Print summary statistics
cat("\n=== Linkage Count Summary (Per Sample, Filtered) ===\n\n")
cat("By Environment (per sample, ≥5 samples):\n")
print(env_sample_counts %>% pivot_wider(names_from = MGE_type, values_from = count, values_fill = 0))

cat("\nBy Phylum (Top 10, per sample, ≥5 samples):\n")
print(phylum_sample_counts %>% pivot_wider(names_from = MGE_type, values_from = count, values_fill = 0))

cat("\nBoxplots saved to linkage_counts_env_boxplot.pdf and linkage_counts_phylum_boxplot.pdf\n")

degree_df <- read_csv("/home/shuaiw/mGlu/tmp/figures/multi_env_linkage/network_99/virus_plasmid_degrees.csv", show_col_types = FALSE)

# Prepare degree data (only plasmid and virus)
degree_df <- degree_df %>%
  filter(type %in% c("plasmid", "virus")) %>%
  mutate(type = factor(type, levels = c("plasmid", "virus")))

# Bar plot: mean ± SD of degree by MGE type
# Compute summary (mean, sd, n)
deg_summary <- degree_df %>%
  group_by(type) %>%
  summarise(mean_deg = mean(degree, na.rm = TRUE),
            sd_deg = sd(degree, na.rm = TRUE),
            n = n(), .groups = "drop")

# Compute t-test p-value for annotation
tt <- try(t.test(degree ~ type, data = degree_df), silent = TRUE)
if (!inherits(tt, "try-error")) {
  pval_text <- format.pval(tt$p.value, digits = 3, eps = 1e-4)
} else {
  pval_text <- "NA"
}

p_deg <- ggplot(deg_summary, aes(x = type, y = mean_deg, fill = type)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean_deg - sd_deg, ymax = mean_deg + sd_deg), width = 0.2) +
  scale_fill_discrete() +
  labs(x = "", y = "Degree", title = "") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)) +
  annotate("text", x = 1.5, y = max(deg_summary$mean_deg + deg_summary$sd_deg, na.rm = TRUE) * 1.05, label = paste0("p = ", pval_text), size = 5) +
  theme(legend.position = "none")

ggsave("../../tmp/figures/multi_env_linkage/network_99/virus_plasmid_degree_barplot.pdf", plot = p_deg, width = 3, height = 5)

# Print t-test result
if (!inherits(tt, "try-error")) {
  cat(sprintf("\nT-test (degree) plasmid vs virus: p = %.4g\n", tt$p.value))
} else {
  cat("\nT-test (degree) failed or insufficient data.\n")
}
