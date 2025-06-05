library(ggplot2)
library(readr)
library(gridExtra)
library(cowplot)
library(tidyverse) 


# Read the CSV file
df <- read_csv("/home/shuaiw/borg/paper/linkage/subsample_96plex.csv")
df2 <- read_csv("/home/shuaiw/borg/paper/linkage/subsample_96plex_depth.csv")
df
df$recall


# Convert proportion to numeric and sort
df <- df %>%
  mutate(proportion = as.numeric(proportion)) %>%
  arrange(proportion)

df2 <- df2 %>%
  mutate(proportion = as.numeric(proportion)) %>%
  arrange(proportion)

## covert depth_cutoff to factor
df$depth_cutoff <- factor(df$depth_cutoff)

# Create the line plot with color representing coverage (depth_cutoff)
p1 <- ggplot(df, aes(x = proportion, y = recall, color = depth_cutoff, group = depth_cutoff)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "Proportion", y = "Recall", color = "Depth Cutoff")

p2 <- ggplot(df, aes(x = proportion, y = precision, color = depth_cutoff, group = depth_cutoff)) +
  geom_line() +
  geom_point(size = 2) +
  theme_bw() +
  labs(x = "Proportion", y = "Precision", color = "Depth Cutoff")




# # Combine plots without legends
prow <- plot_grid(
  p1 + theme(legend.position = "bottom"),
  p2 + theme(legend.position = "none"),
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
pdf(file = "../../tmp/results/96plex_subsample.pdf", width = 7, height = 4, onefile = FALSE)
  print(combined_plot)
dev.off()




