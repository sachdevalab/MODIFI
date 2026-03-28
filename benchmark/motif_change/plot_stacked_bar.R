paper_fig_dir = "../../tmp/figures/strain_diff/drep_99/"
# paper_fig_dir = "../../tmp/figures/strain_diff/iso_drep_99/"
stack_df <- read.csv(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.csv"))
# Sort phylum so that 'others' is at the bottom, and the rest by total count (descending)
phylum_counts <- aggregate(count ~ phylum, data = stack_df, sum)
phylum_counts <- phylum_counts[order(-phylum_counts$count), ]
if ("others" %in% phylum_counts$phylum) {
  phylum_levels <- c(as.character(phylum_counts$phylum[phylum_counts$phylum != "others"]), "others")
} else {
  phylum_levels <- as.character(phylum_counts$phylum)
}
# Reverse for legend so 'others' is at the bottom
stack_df$phylum <- factor(stack_df$phylum, levels = rev(phylum_levels))
# Remove 'p__' prefix from phylum names for legend display
phylum_labels <- gsub('^p__', '', levels(stack_df$phylum))
library(ggplot2)
library(ggpattern)

# Create pattern column: solid for normal, stripe for variation
stack_df$pattern <- ifelse(stack_df$motif_variation_flag == "variation", "stripe", "none")

ggplot(stack_df, aes(x=as.factor(cutoff), y=count, fill=phylum)) +
  geom_bar_pattern(aes(pattern = pattern), 
                   stat="identity", 
                   position="stack",
                   width = 0.7,
                   pattern_density = 0.1,
                   pattern_spacing = 0.02,
                   pattern_angle = 45) +
  scale_fill_manual(
    values = setNames(RColorBrewer::brewer.pal(length(levels(stack_df$phylum)), "Set2"), levels(stack_df$phylum)),
    labels = phylum_labels
  ) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"),
                       labels = c("none" = "Normal", "stripe" = "Variation"),
                       name = "") +
  xlab("Cluster member cutoff") +
  ylab("No. of clusters") +
  labs(fill = "Phylum") +
  guides(pattern = "none") +  ## control whether to show pattern legend
  theme_minimal() +
    theme(axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 10))
ggsave(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.pdf"), width = 5, height = 3)
# ggsave(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.pdf"), width = 5, height = 3)