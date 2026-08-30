#!/usr/bin/env Rscript
# Figure 4a-b (ECE-host motif sharing) for the isolate set, run on BOTH the ORIGINAL full linked set
# (1124) and the NEW strict-filtered set (filterpass_isolate_FINAL, pass==1 -> 870).
#   panel a = proportion of ECEs with an IDENTICAL host motif set (jaccard_similarity_filtered == 1),
#             by phylum x ECE type
#   panel b = mean jaccard_similarity_filtered +/- SE, same grouping
# Reuses the plotting logic of benchmark/ece_anno/plot_identity_veryhigh.R. Nothing existing is overwritten;
# all outputs go to tmp/rev_figs/ece_anno/isolate/.
#   Rscript: /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_identity_strict.R
suppressMessages({library(dplyr); library(ggplot2); library(readr); library(gridExtra)})

JACCARD_CSV  <- "/home/shuaiw/MODIFI/tmp/figures/motif_sharing/jaccard_same_sample.csv"
STRICT_CSV   <- "/home/shuaiw/borg/revision/ece_anno/isolate_all/filterpass_isolate_FINAL.csv"
OUT_DIR      <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/isolate"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

make_fig <- function(df, tag, data_csv, fig_pdf, fig_png) {
  om <- df %>% filter(jaccard_similarity == 1)
  fm <- df %>% filter(jaccard_similarity_filtered == 1)
  cat(sprintf("\n[%s] ECEs = %d  (plasmid %d, virus %d)\n", tag, nrow(df),
              sum(df$mge_type == "plasmid"), sum(df$mge_type == "virus")))
  cat(sprintf("[%s] identical (jaccard==1):          %.4f (%d/%d)\n", tag,
              nrow(om)/nrow(df), nrow(om), nrow(df)))
  cat(sprintf("[%s] identical (jaccard_filtered==1): %.4f (%d/%d)\n", tag,
              nrow(fm)/nrow(df), nrow(fm), nrow(df)))
  write.csv(df, data_csv, row.names = FALSE)

  combos <- expand.grid(phylum = unique(df$phylum), mge_type = unique(df$mge_type),
                        stringsAsFactors = FALSE)

  # panel a: proportion identical (filtered)
  a <- df %>% group_by(phylum, mge_type) %>%
    summarise(perfect_count = sum(jaccard_similarity_filtered == 1), total_count = n(),
              proportion = perfect_count / total_count, .groups = "drop")
  a <- combos %>% left_join(a, by = c("phylum", "mge_type")) %>%
    mutate(across(c(perfect_count, total_count, proportion), ~ ifelse(is.na(.), 0, .)))
  p4 <- ggplot(a, aes(phylum, proportion, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = total_count), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) + labs(x = "Phylum", y = "Proportion of identical motifs", fill = "ECE type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11), axis.title = element_text(size = 12),
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.major.x = element_blank(), legend.position = "top")

  # panel b: mean jaccard (filtered) +/- SE
  bstat <- df %>% group_by(phylum, mge_type) %>%
    summarise(mean = mean(jaccard_similarity_filtered, na.rm = TRUE),
              se = sd(jaccard_similarity_filtered, na.rm = TRUE) / sqrt(n()),
              count = n(), .groups = "drop")
  bstat <- combos %>% left_join(bstat, by = c("phylum", "mge_type")) %>%
    mutate(across(c(mean, se, count), ~ ifelse(is.na(.), 0, .)))
  p6 <- ggplot(bstat, aes(phylum, mean, fill = mge_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5) +
    geom_text(aes(label = count, y = mean + se), position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3, color = "black") +
    ylim(0, 1) + labs(x = "Phylum", y = "Mean Jaccard Similarity", fill = "ECE type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11), axis.title = element_text(size = 12),
          panel.grid.major.y = element_line(colour = "grey80"),
          panel.grid.major.x = element_blank(), legend.position = "top")

  shared <- get_legend(p4)
  combined <- arrangeGrob(
    arrangeGrob(p4 + theme(legend.position = "none"),
                p6 + theme(legend.position = "none"), ncol = 2),
    shared, nrow = 2, heights = c(10, 0.5))
  ggsave(fig_pdf, combined, width = 7, height = 5, dpi = 400)
  ggsave(fig_png, combined, width = 7, height = 5, dpi = 200)
}

full <- read.csv(JACCARD_CSV)
strict_keys <- read_csv(STRICT_CSV, show_col_types = FALSE) %>%
  filter(pass == 1) %>% select(sample, MGE)
strict <- full %>% inner_join(strict_keys, by = c("prefix" = "sample", "mge_contig" = "MGE"))

make_fig(full,   "ORIGINAL full",
         file.path(OUT_DIR, "jaccard_same_sample_full.csv"),
         file.path(OUT_DIR, "jaccard_identity_full_isolate.pdf"),
         file.path(OUT_DIR, "jaccard_identity_full_isolate.png"))
make_fig(strict, "STRICT (pass==1)",
         file.path(OUT_DIR, "jaccard_same_sample_strict.csv"),
         file.path(OUT_DIR, "jaccard_identity_strict_isolate.pdf"),
         file.path(OUT_DIR, "jaccard_identity_strict_isolate.png"))
cat("\nSaved full + strict figures + CSVs to", OUT_DIR, "\n")
