#!/usr/bin/env Rscript
# Figure 4a-b (motif sharing between ECEs and their hosts) restricted to the
# very-high-confidence isolate ECE subset.
#
# Reuses the plotting logic of benchmark/isolation/plot_identity_motifs.R:
#   panel a = proportion of ECEs with an IDENTICAL host motif set
#             (jaccard_similarity_filtered == 1), by phylum x MGE type
#   panel b = mean jaccard_similarity_filtered +/- SE, same grouping
#
# Nothing existing is overwritten. New outputs (both in tmp/rev_figs/ece_anno/):
#   jaccard_same_sample_veryhigh.csv       (filtered data used to make the plot)
#   jaccard_identity_veryhigh_isolate.pdf  (the figure)
#
# Run with the r_env conda R:
#   /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_identity_veryhigh.R

library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)

JACCARD_CSV <- "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"
EVIDENCE_TSV <- "/home/shuaiw/borg/revision/ece_anno/ece_evidence_all.tsv"
OUT_DIR <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"

plot_veryhigh <- function() {
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  # --- load + filter to the very-high-confidence ECE subset ---
  full_df <- read.csv(JACCARD_CSV)
  evid <- read_tsv(EVIDENCE_TSV, show_col_types = FALSE) %>%
    select(sample, seq_name, very_high_confidence)

  same_sample_df <- full_df %>%
    inner_join(evid, by = c("prefix" = "sample", "mge_contig" = "seq_name")) %>%
    filter(very_high_confidence == TRUE)

  cat(sprintf("Very-high-confidence ECEs retained: %d / %d\n",
              nrow(same_sample_df), nrow(full_df)))

  # save the exact table used to make the plot
  write.csv(same_sample_df,
            file.path(OUT_DIR, "jaccard_same_sample_veryhigh.csv"),
            row.names = FALSE)

  # --- overall proportions (identical = jaccard_similarity_filtered == 1) ---
  om <- same_sample_df %>% filter(jaccard_similarity == 1)
  fm <- same_sample_df %>% filter(jaccard_similarity_filtered == 1)
  cat(sprintf("Overall identical (jaccard_similarity==1):          %.4f (%d/%d)\n",
              nrow(om)/nrow(same_sample_df), nrow(om), nrow(same_sample_df)))
  cat(sprintf("Overall identical (jaccard_similarity_filtered==1): %.4f (%d/%d)\n",
              nrow(fm)/nrow(same_sample_df), nrow(fm), nrow(same_sample_df)))

  all_phylums <- unique(same_sample_df$phylum)
  all_mge_types <- unique(same_sample_df$mge_type)
  all_combinations <- expand.grid(phylum = all_phylums, mge_type = all_mge_types,
                                  stringsAsFactors = FALSE)

  # --- panel a: proportion identical (filtered) by phylum x mge_type ---
  phylum_mge_df2 <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(perfect_count = sum(jaccard_similarity_filtered == 1),
              total_count = n(),
              proportion = perfect_count / total_count, .groups = "drop")
  phylum_mge_df2 <- all_combinations %>%
    left_join(phylum_mge_df2, by = c("phylum", "mge_type")) %>%
    mutate(perfect_count = ifelse(is.na(perfect_count), 0, perfect_count),
           total_count = ifelse(is.na(total_count), 0, total_count),
           proportion = ifelse(is.na(proportion), 0, proportion))

  cat("\n--- Proportion identical (filtered) by Phylum x MGE Type ---\n")
  for (i in 1:nrow(phylum_mge_df2)) {
    cat(sprintf("%s - %s: %.4f (%d/%d)\n", phylum_mge_df2$phylum[i],
                phylum_mge_df2$mge_type[i], phylum_mge_df2$proportion[i],
                phylum_mge_df2$perfect_count[i], phylum_mge_df2$total_count[i]))
  }

  p4 <- ggplot(phylum_mge_df2, aes(x = phylum, y = proportion, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = total_count), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(x = "Phylum", y = "Proportion of identical motifs", fill = "ECE type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.major.x = element_blank(),
          legend.position = "top")

  # --- panel b: mean jaccard_similarity_filtered +/- SE ---
  jaccard_filtered_stats <- same_sample_df %>%
    group_by(phylum, mge_type) %>%
    summarise(mean = mean(jaccard_similarity_filtered, na.rm = TRUE),
              se = sd(jaccard_similarity_filtered, na.rm = TRUE) / sqrt(n()),
              count = n(), .groups = "drop")
  jaccard_filtered_stats <- all_combinations %>%
    left_join(jaccard_filtered_stats, by = c("phylum", "mge_type")) %>%
    mutate(mean = ifelse(is.na(mean), 0, mean),
           se = ifelse(is.na(se), 0, se),
           count = ifelse(is.na(count), 0, count))

  cat("\n--- Mean Jaccard Similarity (filtered) by Phylum x MGE Type ---\n")
  for (i in 1:nrow(jaccard_filtered_stats)) {
    cat(sprintf("%s - %s: mean = %.4f +/- %.4f (n=%d)\n",
                jaccard_filtered_stats$phylum[i], jaccard_filtered_stats$mge_type[i],
                jaccard_filtered_stats$mean[i], jaccard_filtered_stats$se[i],
                jaccard_filtered_stats$count[i]))
  }

  p6 <- ggplot(jaccard_filtered_stats, aes(x = phylum, y = mean, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5) +
    geom_text(aes(label = count, y = mean + se), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) +
    labs(x = "Phylum", y = "Mean Jaccard Similarity", fill = "ECE type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.major.x = element_blank(),
          legend.position = "top")

  get_legend <- function(myggplot) {
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    tmp$grobs[[leg]]
  }

  shared_legend <- get_legend(p4)
  combined_plot_46 <- arrangeGrob(
    arrangeGrob(p4 + theme(legend.position = "none"),
                p6 + theme(legend.position = "none"), ncol = 2),
    shared_legend, nrow = 2, heights = c(10, 0.5))

  ggsave(file.path(OUT_DIR, "jaccard_identity_veryhigh_isolate.pdf"),
         combined_plot_46, width = 7, height = 5, dpi = 400)
  ggsave(file.path(OUT_DIR, "jaccard_identity_veryhigh_isolate.png"),
         combined_plot_46, width = 7, height = 5, dpi = 200)
  cat("\nSaved figure (pdf + png) + filtered CSV to", OUT_DIR, "\n")
}

plot_veryhigh()
