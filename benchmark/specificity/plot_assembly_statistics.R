library(readr)
library(ggplot2)
library(patchwork)

df <- read_csv("/home/shuaiw/mGlu/tmp/figures/multi_env_linkage/genome_data_all_samples.csv")

# Define color map for environments
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

# Create N50 box plot
p1 <- ggplot(df, aes(x = environment, y = N50, fill = environment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       x = "Environment",
       y = "N50 (bp)") +
  scale_y_log10() +
  scale_fill_manual(values = env_colors)

# Create genome size box plot
p2 <- ggplot(df, aes(x = environment, y = genome_size, fill = environment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       x = "Environment",
       y = "Genome Size (bp)") +
  scale_y_log10() +
  scale_fill_manual(values = env_colors)

# Combine plots
combined_plot <- p1 + p2

# Display the plot
print(combined_plot)

# Save the plot
ggsave("/home/shuaiw/mGlu/tmp/figures/multi_env_linkage/assembly_statistics_by_environment.pdf", 
       combined_plot, 
       width = 8, 
       height = 4)
