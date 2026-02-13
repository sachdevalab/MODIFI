library(ggplot2)
library(dplyr)

plot_MTase <- function(df_all_data, fig_dir) {
  # Define custom color palette for environments
  env_colors <- c(
    'mock' = '#e41a1c',
    'mice_gut' = '#377eb8',
    'sugarcane' = '#4daf4a',
    'infant_gut' = '#984ea3',
    'cow_rumen' = '#ff7f00',
    'cow_bioreactor' = '#ffff33',
    'adult_gut' = '#a65628',
    'ocean' = '#f781bf',
    'soil' = '#999999'
  )
  
  # Calculate Pearson correlation
  cor_test <- cor.test(df_all_data$RM_num, df_all_data$motif_num, method = "pearson")
  r_value <- cor_test$estimate
  p_value <- cor_test$p.value
  

    p_text <- sprintf("p = %.2e", p_value)

  annotation_text <- sprintf("r = %.3f\n%s", r_value, p_text)
  
  # Scatter plot with regression line and jitter to handle overlapping points
  p1 <- ggplot(df_all_data, aes(x = RM_num, y = motif_num, color = environment)) +
    geom_jitter(alpha = 0.5, size = 2.5, width = 0.2, height = 0.2) +
    geom_smooth(aes(x = RM_num, y = motif_num), 
                method = "lm", se = TRUE, color = "black", 
                fill = "gray80", alpha = 0.3,
                inherit.aes = FALSE) +
    scale_color_manual(values = env_colors) +
    labs(title = "",
         x = "No. of MTases",
         y = "No. of Motifs",
         color = "Environment") +
    annotate("text", x = Inf, y = Inf, label = annotation_text, 
             hjust = 1.1, vjust = 1.5, size = 4.5, fontface = "plain") +
    theme_minimal() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.title = element_text(size = 12))
  
  # Save plot
  output_file <- file.path(fig_dir, "MTase_vs_motif_num.pdf")
  ggsave(filename = output_file, 
         plot = p1, 
         width = 5, height = 4)
  
  cat("Plot saved to:", output_file, "\n")
  cat(sprintf("Pearson correlation: r = %.3f, p = %.3e\n", r_value, p_value))
  
  return(p1)
}

# Example usage:
fig_dir <- "../../tmp/figures/multi_env_linkage/"
df_all_data <- read.csv(file.path(fig_dir, "motif_num_all_samples.csv"))

# Generate plot
plot_MTase(df_all_data, fig_dir)
