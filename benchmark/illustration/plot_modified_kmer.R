library(dplyr)
library(ggplot2)
library(readr)

# Assuming 'data' is a list or matrix with kmer and ipd columns
# df <- data.frame(kmer = data[,1], ipd = data[,2])
# Or if you already have a dataframe, just use it directly


df <- read_csv("/home/shuaiw/borg/paper/workflow_fig/kmer_ipd_distribution.csv")


# Filter data for the specific kmer GCGCGCATCG and contig 96plex_11_C
df_selected <- df %>%
  filter(kmer == "GCCGGCCTGC",
         grepl("96plex_21_C", contig, ignore.case = TRUE))

cat("Number of data points for GCGCGCATCG from 96plex_11_C:", nrow(df_selected), "\n")

# Plot boxplot with jitter points for 96plex_11_C only
p2 <- ggplot(df_selected, aes(x = "", y = ipd)) +
  geom_boxplot(fill = "lightgray", color = "black", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, color = "steelblue") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = 11),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 13, face = "bold"),
        axis.ticks.x = element_blank()) +
  labs(x = "", 
       y = "IPD",
       title = "")

# Save the density plot
ggsave("/home/shuaiw/mGlu/tmp/figures/framework/kmer_ipd_modified.png", 
       plot = p2, 
       width = 3, 
       height = 3, 
       units = "in", dpi = 400)

# Print statistics for the kmer
cat("\nStatistics for GCGCGCATCG (96plex_11_C):\n")
cat("Mean IPD:", mean(df_selected$ipd, na.rm = TRUE), "\n")
cat("SD IPD:", sd(df_selected$ipd, na.rm = TRUE), "\n")
