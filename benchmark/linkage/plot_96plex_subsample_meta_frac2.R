library(ggplot2)
library(readr)
library(gridExtra)
library(cowplot)
library(tidyverse) 


# Read the CSV file
df <- read_csv("/home/shuaiw/borg/paper/linkage/subsample_96plex.csv")
df_dp <- read_csv("/home/shuaiw/borg/paper/linkage/subsample_96plex_meta_depth.csv")
df
df$recall


# Convert fraction to numeric and sort
df <- df %>%
  mutate(fraction = as.numeric(fraction)) %>%
  arrange(fraction)


## covert depth_cutoff to factor
df$depth_cutoff <- factor(df$depth_cutoff)

# Create the line plot with color representing coverage (depth_cutoff)
p1 <- ggplot(df, aes(x = fraction, y = recall)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "96plex fraction (%)", y = "Recall")

p2 <- ggplot(df, aes(x = fraction, y = precision)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "96plex fraction (%)", y = "Precision")

## plot boxplot for depth, each fraction has a boxplot
df_dp <- df_dp %>%
  mutate(fraction = as.numeric(fraction)) %>%
  arrange(fraction)
df_dp
df_dp$fraction <- factor(df_dp$fraction)
p3 <- ggplot(df_dp, aes(x = fraction, y = depth)) +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  labs(x = "96plex fraction (%)", y = "Depth") +
  ylim(0, 150)


# # Combine plots without legends
prow <- plot_grid(
  p1 + theme(legend.position = "bottom"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

# Extract the legend from one of the plots (without removing it)
legend <- get_legend(
  p1 +
    guides(color = guide_legend(nrow = 1)) +  # Ensure the legend is properly formatted
    theme(legend.position = "bottom")
)

# Combine the plots and the legend
combined_plot <- plot_grid(prow, legend, ncol = 1, rel_heights = c(1, 0.2))

# Now write to PDF
pdf(file = "../../tmp/results/96plex_subsample_frac.pdf", width = 10, height = 4, onefile = FALSE)
  print(combined_plot)
dev.off()




