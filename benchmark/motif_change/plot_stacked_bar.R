# paper_fig_dir = "../../tmp/figures/strain_diff/drep_99/"
paper_fig_dir = "../../tmp/figures/strain_diff/iso_drep_99/"
stack_df <- read.csv(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.csv"))
library(ggplot2)
library(ggpattern)

# Create pattern column: solid for normal, stripe for variation
stack_df$pattern <- ifelse(stack_df$motif_variation_flag == "variation", "stripe", "none")

ggplot(stack_df, aes(x=as.factor(cutoff), y=count, fill=phylum)) +
  geom_bar_pattern(aes(pattern = pattern), 
                   stat="identity", 
                   position="stack",
                   pattern_density = 0.1,
                   pattern_spacing = 0.02,
                   pattern_angle = 45) +
  scale_fill_brewer(palette = "Set2") +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"),
                       labels = c("none" = "Normal", "stripe" = "Variation"),
                       name = "Variation Flag") +
  xlab("Strain member cutoff") +
  ylab("No. of strains") +
  labs(fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "right")
ggsave(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.pdf"), width = 5, height = 3)
# ggsave(paste0(paper_fig_dir, "motif_variation_counts_by_phylum.pdf"), width = 5, height = 3)