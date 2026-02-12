library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)

df <- read_csv("/home/shuaiw/borg/paper/run2/96plex/96plex_methylation4/ipd_ratio/96plex_1_C.ipd3.csv")

# Calculate mean and standard deviation for ipd_ratio
mean_ipd <- mean(df$ipd_ratio, na.rm = TRUE)
sd_ipd <- sd(df$ipd_ratio, na.rm = TRUE)

# Calculate threshold for p-value = 0.001 (right tail, 99th percentile)
threshold <- mean_ipd + qnorm(0.999) * sd_ipd

# Count significantly higher points
n_significant <- sum(df$ipd_ratio > threshold, na.rm = TRUE)

cat("IPD Ratio Statistics:\n")
cat("Mean:", mean_ipd, "\n")
cat("SD:", sd_ipd, "\n")
cat("Threshold (p < 0.001, right tail):", threshold, "\n")
cat("Number of significantly higher points:", n_significant, "\n")
cat("Percentage:", round(n_significant / nrow(df) * 100, 2), "%\n\n")

# Plot distribution of tMean
p1 <- ggplot(df, aes(x = tMean)) +
  geom_density(fill = "steelblue", color = "steelblue", alpha = 0.4, size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  labs(x = "IPD", 
       y = "Density",
       title = "")

# Plot distribution of ipd_ratio with threshold line
p2 <- ggplot(df, aes(x = ipd_ratio)) +
  geom_density(fill = "coral", color = "coral", alpha = 0.4, size = 1) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "black", size = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11)) +
  labs(x = "IPD Ratio", 
       y = "Density",
       title = "")

# Combine the two plots
combined_plot <- grid.arrange(p1, p2, ncol = 2)

# Save the combined plot
ggsave("/home/shuaiw/mGlu/tmp/figures/framework/ipd_ratio.png", 
       plot = combined_plot, 
       width = 6, 
       height = 3, 
       units = "in", 
       dpi = 400)

# Print summary statistics
cat("Summary statistics for tMean:\n")
print(summary(df$tMean))

cat("\nSummary statistics for ipd_ratio:\n")
print(summary(df$ipd_ratio))
