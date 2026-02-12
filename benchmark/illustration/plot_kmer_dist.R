library(dplyr)
library(ggplot2)
library(readr)

# Assuming 'data' is a list or matrix with kmer and ipd columns
# df <- data.frame(kmer = data[,1], ipd = data[,2])
# Or if you already have a dataframe, just use it directly


df <- read_csv("/home/shuaiw/borg/paper/workflow_fig/kmer_ipd_distribution.csv")


# Select 4 kmers with the most different mean IPD values
# First calculate mean and std for each kmer
kmer_stats <- df %>%
  group_by(kmer) %>%
  summarise(mean_ipd = mean(ipd, na.rm = TRUE),
            sd_ipd = sd(ipd, na.rm = TRUE)) %>%
  arrange(mean_ipd)

# Remove kmers with too high standard deviation (top 25% by std)
sd_threshold <- quantile(kmer_stats$sd_ipd, 0.75, na.rm = TRUE)
kmer_means <- kmer_stats %>%
  filter(sd_ipd <= sd_threshold)

## echo number of kmers
cat("Number of kmers after filtering by std:", nrow(kmer_means), "\n")
# Select kmers with maximum spread: lowest, 33rd percentile, 67th percentile, and highest
n_kmers <- nrow(kmer_means)
if (n_kmers >= 4) {
  selected_indices <- c(1, 
                        round(n_kmers /7), 
                        round(n_kmers * 0.67), 
                        round(n_kmers * 1))
  selected_kmers <- kmer_means$kmer[selected_indices]
} else {
  selected_kmers <- kmer_means$kmer
}

# Filter data for selected kmers
df_selected <- df %>%
  filter(kmer %in% selected_kmers) %>%
  mutate(kmer = factor(kmer, levels = selected_kmers))

# Plot density distributions in a single plot with different colors
p2 <- ggplot(df_selected, aes(x = ipd, fill = kmer, color = kmer)) +
  geom_density(alpha = 0.4, size = 1) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.65, 0.55),
        legend.title = element_text(size = 11 ),
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  labs(x = "IPD", 
       y = "Density",
       fill = "K-mer",
       color = "K-mer",
       title = "")

# Save the density plot
ggsave("/home/shuaiw/mGlu/tmp/figures/framework/kmer_ipd_density_selected.png", 
       plot = p2, 
       width = 4, 
       height = 3, 
       units = "in", dpi = 400)

# Print selected kmers and their mean IPD values
print("Selected kmers and their mean IPD values:")
print(kmer_means %>% filter(kmer %in% selected_kmers))
