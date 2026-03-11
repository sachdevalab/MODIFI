#!/usr/bin/env Rscript

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

# Read data
gc_df <- read_csv("../../tmp/figures/multi_env_linkage/network_99/mge_host_gc_cov.csv", show_col_types = FALSE)

# Load phylum -> Gram-Positive / Gram-Negative / Archaea (canonical: benchmark/specificity/phylum_gram_archaea_annotation.csv)
annotation_file <- NULL
for (path in c("../specificity/phylum_gram_archaea_annotation.csv",
               "benchmark/specificity/phylum_gram_archaea_annotation.csv",
               "phylum_gram_archaea_annotation.csv",
               "../../benchmark/specificity/phylum_gram_archaea_annotation.csv")) {
  if (file.exists(path)) {
    annotation_file <- path
    break
  }
}
if (is.null(annotation_file)) stop("Cannot find phylum_gram_archaea_annotation.csv")
phylum_anno <- read_csv(annotation_file, show_col_types = FALSE)
phylum_class_map <- setNames(phylum_anno$classification, phylum_anno$phylum)

gc_df <- gc_df %>%
  mutate(MGE_type = factor(MGE_type, levels = c("plasmid", "virus"))) %>%
  mutate(host_phylum_clean = gsub("^p__", "", host_phylum)) %>%
  mutate(domain_type = ifelse(host_phylum_clean %in% names(phylum_class_map),
                              phylum_class_map[host_phylum_clean],
                              "Gram-Negative")) %>%
  mutate(domain_type = factor(domain_type, levels = c("Gram-Positive", "Gram-Negative", "Archaea")))

## count the mean cos_sim  and std for all linkages, and for plasmid and virus separately
gc_summary <- gc_df %>%
  group_by(MGE_type) %>%
  summarise(mean_cos_sim = mean(cos_sim, na.rm = TRUE),
            sd_cos_sim = sd(cos_sim, na.rm = TRUE),
            n = n(), .groups = "drop")
cat("=== Cosine Similarity Summary ===\n")
print(gc_summary)
## also print the overall mean and sd
overall_summary <- gc_df %>%
  summarise(mean_cos_sim = mean(cos_sim, na.rm = TRUE),
            sd_cos_sim = sd(cos_sim, na.rm = TRUE),
            n = n())
cat("\nOverall Cosine Similarity:\n")
print(overall_summary)

# Count linkages per sample by environment and MGE type
env_sample_counts <- gc_df %>%
  group_by(environment, sample, MGE_type, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")

# Only retain environments with >=2 samples (for boxplot; cow_bioreactor etc. included)
env_sample_counts <- env_sample_counts %>%
  group_by(environment) %>%
  filter(n_distinct(sample) >= 4) %>%
  ungroup()

# Order environments by total count (for consistent x-axis order)
env_order <- gc_df %>%
  count(environment, sort = TRUE) %>%
  pull(environment)
env_sample_counts$environment <- factor(env_sample_counts$environment, levels = env_order)

# Count linkages per sample by phylum and MGE type (top 10 phyla with >=5 samples)
# First filter phyla with >=5 samples, then select top 10 by count
phyla_with_enough_samples <- gc_df %>%
  group_by(host_phylum) %>%
  filter(n_distinct(sample) >= 5) %>%
  ungroup()

top_phyla <- phyla_with_enough_samples %>%
  count(host_phylum, sort = TRUE) %>%
  head(7) %>%
  pull(host_phylum)

phylum_sample_counts <- phyla_with_enough_samples %>%
  filter(host_phylum %in% top_phyla) %>%
  mutate(host_phylum = gsub("^p__", "", host_phylum)) %>%
  group_by(host_phylum, sample, MGE_type, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")
phylum_sample_counts$host_phylum <- factor(phylum_sample_counts$host_phylum, levels = gsub("^p__", "", top_phyla))

# Count linkages per sample by domain type (Gram-Positive / Gram-Negative / Archaea) and MGE type
domain_sample_counts <- gc_df %>%
  group_by(domain_type, sample, MGE_type, .drop = FALSE) %>%
  summarise(count = n(), .groups = "drop")

# Optional: normalize by MGE count per sample (from count_mge_per_sample.py; same logic as profile_good_ctgs.py)
paper_fig_dir <- "../../tmp/figures/multi_env_linkage/network_99"
mge_counts_file <- file.path(paper_fig_dir, "mge_counts_per_sample.csv")
domain_sample_counts_norm <- NULL
if (file.exists(mge_counts_file)) {
  mge_counts <- read_csv(mge_counts_file, show_col_types = FALSE)
  domain_sample_counts_norm <- domain_sample_counts %>%
    left_join(mge_counts %>% select(sample, n_plasmid, n_virus), by = "sample") %>%
    mutate(
      n_mge = if_else(MGE_type == "plasmid", n_plasmid, n_virus),
      n_mge = replace(n_mge, !is.na(n_mge) & n_mge == 0, 1),
      count_norm = count / n_mge
    ) %>%
    filter(!is.na(count_norm))
} else {
  cat("Note: mge_counts_per_sample.csv not found. Run: python count_mge_per_sample.py\n")
}

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
    x = "Habitat",
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

# Boxplot 3: Linkages per sample by domain type (Gram-Positive / Gram-Negative / Archaea), hue = MGE type
domain_pvals <- domain_sample_counts %>%
  group_by(domain_type) %>%
  filter(length(unique(MGE_type)) == 2) %>%
  summarise(p = tryCatch(t.test(count ~ MGE_type)$p.value, error = function(e) NA_real_)) %>%
  mutate(star = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "")))
domain_ypos <- domain_sample_counts %>% group_by(domain_type) %>% summarise(y = max(count, na.rm = TRUE) * 1.1)
domain_pvals <- left_join(domain_pvals, domain_ypos, by = "domain_type")

p3 <- ggplot(domain_sample_counts, aes(x = domain_type, y = count, fill = MGE_type)) +
  geom_boxplot(aes(fill = MGE_type), outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 1, color = "black", size = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 1.2, alpha = 1, show.legend = FALSE) +
  scale_fill_discrete(name = "MGE Type") +
  scale_color_discrete(guide = "none") +
  labs(x = "Host classification", y = "Linkages per Sample", title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  geom_text(data = domain_pvals, aes(x = domain_type, y = y, label = star), inherit.aes = FALSE, vjust = 0, size = 6)

ggsave("../../tmp/figures/multi_env_linkage/network_99/linkage_counts_domain_boxplot.pdf", 
       plot = p3, width = 4, height = 5)

# Boxplot 3 (normalized): Linkages per sample per MGE by domain type (requires mge_counts_per_sample.csv)
if (!is.null(domain_sample_counts_norm)) {
  domain_pvals_norm <- domain_sample_counts_norm %>%
    group_by(domain_type) %>%
    filter(length(unique(MGE_type)) == 2) %>%
    summarise(p = tryCatch(t.test(count_norm ~ MGE_type)$p.value, error = function(e) NA_real_)) %>%
    mutate(star = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "")))
  domain_ypos_norm <- domain_sample_counts_norm %>% group_by(domain_type) %>% summarise(y = max(count_norm, na.rm = TRUE) * 1.1)
  domain_pvals_norm <- left_join(domain_pvals_norm, domain_ypos_norm, by = "domain_type")

  p3_norm <- ggplot(domain_sample_counts_norm, aes(x = domain_type, y = count_norm, fill = MGE_type)) +
    geom_boxplot(aes(fill = MGE_type), outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 1, color = "black", size = 0.5) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 1.2, alpha = 1, show.legend = FALSE) +
    scale_fill_discrete(name = "MGE Type") +
    scale_color_discrete(guide = "none") +
    labs(x = "Host classification", y = "Linkages per Sample per MGE", title = "") +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    ) +
    geom_text(data = domain_pvals_norm, aes(x = domain_type, y = y, label = star), inherit.aes = FALSE, vjust = 0, size = 6)

  ggsave(file.path(paper_fig_dir, "linkage_counts_domain_boxplot_norm.pdf"),
         plot = p3_norm, width = 4, height = 5)
  cat("Saved normalized domain boxplot to linkage_counts_domain_boxplot_norm.pdf\n")
}

# Boxplot 4: x = MGE type (plasmid, virus), hue = host type (Gram-Positive vs Gram-Negative; Archaea discarded), normalized, with significance
domain_sample_bacteria <- domain_sample_counts %>%
  filter(domain_type %in% c("Gram-Positive", "Gram-Negative"))
# Use normalized data when available
domain_sample_bacteria_plot <- if (!is.null(domain_sample_counts_norm)) {
  domain_sample_counts_norm %>%
    filter(domain_type %in% c("Gram-Positive", "Gram-Negative"))
} else {
  domain_sample_bacteria %>% mutate(count_norm = count)
}
mge_host_pvals <- domain_sample_bacteria_plot %>%
  group_by(MGE_type) %>%
  filter(n_distinct(domain_type) == 2) %>%
  summarise(p = tryCatch(t.test(count_norm ~ domain_type)$p.value, error = function(e) NA_real_)) %>%
  mutate(star = symnum(p, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", "")))
mge_host_ypos <- domain_sample_bacteria_plot %>% group_by(MGE_type) %>% summarise(y = max(count_norm, na.rm = TRUE) * 1.1)
mge_host_pvals <- left_join(mge_host_pvals, mge_host_ypos, by = "MGE_type")

p4 <- ggplot(domain_sample_bacteria_plot, aes(x = MGE_type, y = count_norm, fill = domain_type)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 1, color = "black", size = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), color = "black", size = 1.2, alpha = 1, show.legend = FALSE) +
  scale_fill_discrete(name = "Host type") +
  labs(x = "MGE Type", y = "Linkages per Sample per MGE", title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  ) +
  geom_text(data = mge_host_pvals, aes(x = MGE_type, y = y, label = star), inherit.aes = FALSE, vjust = 0, size = 6)

ggsave(file.path(paper_fig_dir, "linkage_counts_mge_by_host_type_boxplot.pdf"),
       plot = p4, width = 4, height = 5)

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

# T-tests for domain type
cat("\n=== T-tests: Plasmid vs Virus by Domain (Gram-Positive / Gram-Negative / Archaea) ===\n")
for (dom in levels(domain_sample_counts$domain_type)) {
  dat <- domain_sample_counts %>% filter(domain_type == dom)
  if (length(unique(dat$MGE_type)) == 2) {
    t_res <- try(t.test(count ~ MGE_type, data = dat), silent = TRUE)
    if (!inherits(t_res, "try-error")) {
      cat(sprintf("%s: p = %.4g\n", as.character(dom), t_res$p.value))
    }
  }
}

# T-tests for MGE type panel: Gram-Positive vs Gram-Negative within plasmid and virus (Archaea discarded; normalized when available)
cat("\n=== T-tests: Gram-Positive vs Gram-Negative by MGE type (Archaea discarded) ===\n")
for (mge in c("plasmid", "virus")) {
  dat <- domain_sample_bacteria_plot %>% filter(MGE_type == mge)
  if (n_distinct(dat$domain_type) == 2) {
    t_res <- try(t.test(count_norm ~ domain_type, data = dat), silent = TRUE)
    if (!inherits(t_res, "try-error")) {
      cat(sprintf("%s: p = %.4g\n", mge, t_res$p.value))
    }
  }
}

# Print summary statistics
cat("\n=== Linkage Count Summary (Per Sample, Filtered) ===\n\n")
cat("By Environment (per sample, ≥5 samples):\n")
print(env_sample_counts %>% pivot_wider(names_from = MGE_type, values_from = count, values_fill = 0))

cat("\nBy Phylum (Top 10, per sample, ≥5 samples):\n")
print(phylum_sample_counts %>% pivot_wider(names_from = MGE_type, values_from = count, values_fill = 0))

cat("\nBy Domain (Gram-Positive / Gram-Negative / Archaea, per sample):\n")
print(domain_sample_counts %>% pivot_wider(names_from = MGE_type, values_from = count, values_fill = 0))
if (!is.null(domain_sample_counts_norm)) {
  cat("\nBy Domain (normalized: linkages per sample per MGE):\n")
  print(domain_sample_counts_norm %>% pivot_wider(names_from = MGE_type, values_from = count_norm, values_fill = NA))
}

cat("\nBoxplots saved to linkage_counts_env_boxplot.pdf, linkage_counts_phylum_boxplot.pdf, linkage_counts_domain_boxplot.pdf")
if (!is.null(domain_sample_counts_norm)) cat(", linkage_counts_domain_boxplot_norm.pdf")
cat(", linkage_counts_mge_by_host_type_boxplot.pdf\n")

degree_df <- read_csv("/home/shuaiw/MODIFI/tmp/figures/multi_env_linkage/network_99/virus_plasmid_degrees.csv", show_col_types = FALSE)

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



