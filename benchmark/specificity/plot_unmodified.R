fig_dir <- "../../tmp/figures/multi_env_linkage/"
library(ggplot2)
library(patchwork)
df_all_data <- read.csv(paste0(fig_dir, "motif_num_all_samples.csv"))

# Count genomes with 0 motifs per phylum
zero_motif_genomes <- df_all_data[df_all_data$motif_num == 0, ]
zero_counts <- table(zero_motif_genomes$phylum)
zero_counts_sorted <- sort(zero_counts, decreasing = TRUE)

# Calculate total genomes per phylum
total_counts <- table(df_all_data$phylum)

# Calculate proportion of genomes with 0 motifs per phylum
proportion_df <- data.frame(
  phylum = names(total_counts),
  total_genomes = as.numeric(total_counts),
  zero_motif_genomes = as.numeric(zero_counts[names(total_counts)]),
  stringsAsFactors = FALSE
)
proportion_df$zero_motif_genomes[is.na(proportion_df$zero_motif_genomes)] <- 0
proportion_df$proportion <- proportion_df$zero_motif_genomes / proportion_df$total_genomes
proportion_df$percentage <- proportion_df$proportion * 100

# Sort by proportion (descending)
proportion_df <- proportion_df[order(-proportion_df$proportion), ]

cat("\nProportion of genomes with 0 motifs per phylum:\n")
cat(sprintf("%-40s %10s %15s %12s\n", "Phylum", "Total", "Zero Motifs", "Percentage"))
cat(paste(rep("-", 80), collapse = ""), "\n")
for (i in 1:nrow(proportion_df)) {
  if (proportion_df$proportion[i] > 0) {
    cat(sprintf("%-40s %10d %15d %11.2f%%\n", 
                proportion_df$phylum[i], 
                proportion_df$total_genomes[i], 
                proportion_df$zero_motif_genomes[i], 
                proportion_df$percentage[i]))
  }
}
## print total number of genomes with 0 motifs
cat("\nTotal number of genomes with 0 motifs:", sum(proportion_df$zero_motif_genomes), "\n")
## print total number of genomes 
cat("Total number of genomes:", sum(proportion_df$total_genomes), "\n")

cat("\nNumber of genomes with 0 motifs per phylum:\n")
print(zero_counts_sorted)
cat("\nTotal genomes with 0 motifs:", nrow(zero_motif_genomes), "\n")

# Output RM_num for each genome with 0 motifs
# cat("\n\nRM_num for each genome with 0 motifs:\n")
zero_motif_output <- zero_motif_genomes[, c("contig", "phylum", "RM_num")]
# print(zero_motif_output, row.names = FALSE)

# Summary statistics of RM_num for genomes with 0 motifs
cat("\n\nRM_num summary for genomes with 0 motifs:\n")
cat("Mean:", mean(zero_motif_genomes$RM_num), "\n")
cat("Median:", median(zero_motif_genomes$RM_num), "\n")
cat("Min:", min(zero_motif_genomes$RM_num), "\n")
cat("Max:", max(zero_motif_genomes$RM_num), "\n")
cat("Genomes with RM_num = 0:", sum(zero_motif_genomes$RM_num == 0), "\n")
cat("Genomes with RM_num > 0:", sum(zero_motif_genomes$RM_num > 0), "\n")

# Fisher's exact test for Gram-positive vs Gram-negative enrichment
cat("\n\n" , paste(rep("=", 80), collapse = ""), "\n")
cat("FISHER'S EXACT TEST: Gram-Positive vs Gram-Negative Enrichment\n")
cat("(Only considering phyla with at least one genome without motifs)\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Get phyla that have at least one genome with 0 motifs
phyla_with_zero <- names(zero_counts[zero_counts > 0])

# Filter to only include genomes from these phyla
df_filtered <- df_all_data[df_all_data$phylum %in% phyla_with_zero, ]

cat(sprintf("Total genomes in analysis: %d (from %d phyla with no-motif genomes)\n", 
            nrow(df_filtered), length(phyla_with_zero)))
cat(sprintf("Total genomes excluded: %d (from phyla with no no-motif genomes)\n\n", 
            nrow(df_all_data) - nrow(df_filtered)))

# Clean phylum names and classify by domain/Gram stain
df_filtered$phylum_clean <- gsub("^p__", "", df_filtered$phylum)
gram_positive_phyla <- c("Bacillota", "Bacillota_A", "Bacillota_C", "Bacillota_I", "Actinomycetota")
archaea_phyla <- c("Aenigmatarchaeota", "Methanobacteriota", "Thermoproteota", 
                   "Halobacteriota", "Nanoarchaeota", "Thermoplasmatota")

# Percentage of domain/type among unmodified genomes (motif_num == 0)
zero_motif_genomes$phylum_clean <- gsub("^p__", "", zero_motif_genomes$phylum)
zero_motif_genomes$domain_type <- ifelse(zero_motif_genomes$phylum_clean %in% archaea_phyla, "Archaea",
                                         ifelse(zero_motif_genomes$phylum_clean %in% gram_positive_phyla,
                                                "Gram-Positive", "Gram-Negative"))

zero_total <- nrow(zero_motif_genomes)
zero_gp <- sum(zero_motif_genomes$domain_type == "Gram-Positive")
zero_gn <- sum(zero_motif_genomes$domain_type == "Gram-Negative")
zero_archaea <- sum(zero_motif_genomes$domain_type == "Archaea")

cat("Percentage in unmodified genomes (motif_num == 0):\n")
cat(sprintf("  Gram-Positive: %d / %d (%.2f%%)\n", zero_gp, zero_total, zero_gp / zero_total * 100))
cat(sprintf("  Gram-Negative: %d / %d (%.2f%%)\n", zero_gn, zero_total, zero_gn / zero_total * 100))
cat(sprintf("  Archaea: %d / %d (%.2f%%)\n\n", zero_archaea, zero_total, zero_archaea / zero_total * 100))

df_filtered$domain_type <- ifelse(df_filtered$phylum_clean %in% archaea_phyla, "Archaea",
                                   ifelse(df_filtered$phylum_clean %in% gram_positive_phyla, 
                                          "Gram-Positive", "Gram-Negative"))

# Separate Fisher test: Bacteria only (Gram-Positive vs Gram-Negative)
df_bacteria <- df_filtered[df_filtered$domain_type != "Archaea", ]

cat(sprintf("Bacteria genomes in analysis: %d\n", nrow(df_bacteria)))
cat(sprintf("Archaea genomes in analysis: %d\n\n", sum(df_filtered$domain_type == "Archaea")))

# Create contingency table for bacteria only
# Rows: Gram-Positive, Gram-Negative
# Cols: Zero motifs, Has motifs
gram_pos_zero <- sum(df_bacteria$domain_type == "Gram-Positive" & df_bacteria$motif_num == 0)
gram_pos_with <- sum(df_bacteria$domain_type == "Gram-Positive" & df_bacteria$motif_num > 0)
gram_neg_zero <- sum(df_bacteria$domain_type == "Gram-Negative" & df_bacteria$motif_num == 0)
gram_neg_with <- sum(df_bacteria$domain_type == "Gram-Negative" & df_bacteria$motif_num > 0)

contingency_table <- matrix(c(gram_pos_zero, gram_pos_with, 
                              gram_neg_zero, gram_neg_with),
                           nrow = 2, byrow = TRUE)
rownames(contingency_table) <- c("Gram-Positive", "Gram-Negative")
colnames(contingency_table) <- c("Zero Motifs", "Has Motifs")

cat("Contingency Table:\n")
print(contingency_table)
cat("\n")

# Perform Fisher's exact test
fisher_result <- fisher.test(contingency_table)

cat("Fisher's Exact Test Results:\n")
cat(sprintf("  Odds Ratio: %.4f\n", fisher_result$estimate))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n", fisher_result$conf.int[1], fisher_result$conf.int[2]))
cat(sprintf("  P-value: %.4e\n", fisher_result$p.value))
cat("\n")

# Calculate proportions for interpretation
gram_pos_prop <- gram_pos_zero / (gram_pos_zero + gram_pos_with) * 100
gram_neg_prop <- gram_neg_zero / (gram_neg_zero + gram_neg_with) * 100

cat("Summary:\n")
cat(sprintf("  Gram-Positive genomes without motifs: %d / %d (%.2f%%)\n", 
            gram_pos_zero, gram_pos_zero + gram_pos_with, gram_pos_prop))
cat(sprintf("  Gram-Negative genomes without motifs: %d / %d (%.2f%%)\n", 
            gram_neg_zero, gram_neg_zero + gram_neg_with, gram_neg_prop))
cat("\n")

if (fisher_result$p.value < 0.05) {
  if (fisher_result$estimate > 1) {
    cat("Interpretation: Genomes without motifs are significantly ENRICHED in Gram-Positive bacteria (p < 0.05)\n")
  } else {
    cat("Interpretation: Genomes without motifs are significantly DEPLETED in Gram-Positive bacteria (p < 0.05)\n")
  }
} else {
  cat("Interpretation: No significant difference between Gram-Positive and Gram-Negative bacteria (p >= 0.05)\n")
}
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Publication-style enrichment bar plot (Gram-Positive vs Gram-Negative)
bar_df <- data.frame(
  group = c("Gram-Positive", "Gram-Negative"),
  zero = c(gram_pos_zero, gram_neg_zero),
  total = c(gram_pos_zero + gram_pos_with, gram_neg_zero + gram_neg_with),
  percentage = c(gram_pos_prop, gram_neg_prop),
  stringsAsFactors = FALSE
)
bar_df$group <- factor(bar_df$group, levels = c("Gram-Positive", "Gram-Negative"))
bar_df$label <- sprintf("%d/%d (%.2f%%)", bar_df$zero, bar_df$total, bar_df$percentage)

star_label <- ifelse(fisher_result$p.value < 0.001, "***",
                     ifelse(fisher_result$p.value < 0.01, "**",
                            ifelse(fisher_result$p.value < 0.05, "*", "ns")))
p_label <- paste0("p = ", format(fisher_result$p.value, scientific = TRUE, digits = 2))

y_bracket <- max(23, max(bar_df$percentage) + 3)
y_stars <- y_bracket + 0.45
y_limit <- max(25, y_bracket + 1.0)

p_enrichment <- ggplot(bar_df, aes(x = group, y = percentage, fill = group)) +
  geom_col(width = 0.7, color = "gray25", linewidth = 0.35) +
  annotate("segment", x = 1, xend = 1, y = y_bracket - 1, yend = y_bracket) +
  annotate("segment", x = 1, xend = 2, y = y_bracket, yend = y_bracket) +
  annotate("segment", x = 2, xend = 2, y = y_bracket - 1, yend = y_bracket) +
  annotate("text", x = 1.5, y = y_stars, label = star_label, size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Gram-Positive" = "#808080", "Gram-Negative" = "#808080")) +
  scale_y_continuous(limits = c(0, y_limit), breaks = seq(0, 25, 5), expand = c(0, 0)) +
  labs(
    x = NULL,
    y = "Genomes without motifs (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.grid.major.x = element_blank()
  )

ggsave(paste0(fig_dir, "gram_enrichment_zero_motif_barplot.pdf"),
       p_enrichment, width = 2, height = 4, dpi = 300)
ggsave(paste0(fig_dir, "gram_enrichment_zero_motif_barplot.png"),
       p_enrichment, width = 2, height = 4, dpi = 300)

# Plot bar chart of proportion with zero motifs
# Filter for phyla with proportion > 0
plot_df <- proportion_df[proportion_df$proportion > 0, ]

# Remove p__ prefix from phylum names
plot_df$phylum <- gsub("^p__", "", plot_df$phylum)

# Classify by domain and Gram stain
gram_positive <- c("Bacillota", "Bacillota_A", "Bacillota_C", "Bacillota_I", "Actinomycetota")
archaea <- c("Aenigmatarchaeota", "Methanobacteriota", "Thermoproteota", 
             "Halobacteriota", "Nanoarchaeota", "Thermoplasmatota")
plot_df$domain_type <- ifelse(plot_df$phylum %in% archaea, "Archaea",
                               ifelse(plot_df$phylum %in% gram_positive, 
                                      "Gram-Positive", "Gram-Negative"))

# Sort by number of zero-motif genomes (descending) and set factor levels
plot_df <- plot_df[order(-plot_df$zero_motif_genomes), ]
plot_df$phylum <- factor(plot_df$phylum, levels = plot_df$phylum)

# Plot 1: Number of genomes with 0 motifs
p1 <- ggplot(plot_df, aes(x = phylum, y = zero_motif_genomes, fill = domain_type)) +
  geom_bar(stat = "identity", width = 0.75, alpha = 0.9) +
  geom_text(aes(label = zero_motif_genomes), 
            vjust = -0.5, size = 3.8, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = c("Gram-Positive" = "#8B7AB8", 
                               "Gram-Negative" = "#E85D75",
                               "Archaea" = "#F0A830"),
                    name = "Domain/Type") +
  labs(x = NULL, 
       y = "Number of Genomes\nwithout Motifs") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 11, color = "gray20"),
        axis.title.y = element_text(size = 12, face = "bold", color = "gray10", 
                                     margin = margin(r = 10)),
        legend.position = "none",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.6, "cm"),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 10, 0),
        panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(10, 10, 20, 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)))

# Plot 2: Proportion of genomes with 0 motifs (upside down)
p2 <- ggplot(plot_df, aes(x = phylum, y = -percentage, fill = domain_type)) +
  geom_bar(stat = "identity", width = 0.75, alpha = 0.9) +
  geom_text(data = subset(plot_df, percentage >= 5),
            aes(label = sprintf("%.1f%%", percentage)), 
            vjust = 1.8, size = 3.5, fontface = "bold", color = "gray20") +
  scale_fill_manual(values = c("Gram-Positive" = "#8B7AB8", 
                               "Gram-Negative" = "#E85D75",
                               "Archaea" = "#F0A830"),
                    name = "Domain/Type") +
  labs(x = "Phylum", 
       y = "Percentage without\nMotifs (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, 
                                   color = "gray20", face = "plain"),
        axis.text.y = element_text(size = 11, color = "gray20"),
        axis.title.x = element_text(size = 12, face = "bold", color = "gray10",
                                     margin = margin(t = 10)),
        axis.title.y = element_text(size = 12, face = "bold", color = "gray10",
                                     margin = margin(r = 10)),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(20, 10, 10, 10),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0)),
                     labels = function(x) abs(x))

# Combine plots
combined_plot <- p1 / p2 + 
  plot_layout(heights = c(1.1, 1), guides = "collect") &
  theme(plot.background = element_rect(fill = "white", color = NA),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        legend.key.size = unit(0.6, "cm"),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(0, 0, 10, 0))

ggsave(paste0(fig_dir, "zero_motif_proportion_by_phylum.pdf"), 
       combined_plot, width = 12, height = 9, dpi = 300)
ggsave(paste0(fig_dir, "zero_motif_proportion_by_phylum.png"), 
       combined_plot, width = 12, height = 9, dpi = 300)

print(combined_plot)

print(combined_plot)