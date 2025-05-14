library(ggplot2)
library(readr)
library(gridExtra)
library(cowplot)
library(tidyverse) 


pdf(file="../../tmp/results/96plex_subsample.pdf", width=9, height=4, onefile=FALSE)

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

# Create the line plot with color representing coverage
p1 <- ggplot(df, aes(x = proportion, y = recall)) +
  geom_line() +
  geom_point(size = 2) +
  theme_minimal() +
  labs(x = "Proportion", y = "Recall")

p2 <- ggplot(df, aes(x = proportion, y = precision)) +
  geom_line() +
  ## add point
  geom_point(size = 2) +
  theme_minimal() +
  labs(x = "Proportion", y = "Precision")

df2$proportion <- as.factor(df2$proportion)
## add a boxplot for df2, x=proportion, y=depth
p3 <- ggplot(df2, aes(x = proportion, y = depth )) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Proportion", y = "Depth") +
  ylim(0, 100) 
## print mean and median depth for each proportion
df2 %>%
  group_by(proportion) %>%
  summarise(mean_depth = mean(depth), median_depth = median(depth)) %>%
  print()



# Combine plots without legends
prow <- plot_grid(
  p3 + theme(legend.position = "none"),
  p1 + theme(legend.position = "none"),
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

# Explicitly print the combined plot to ensure it's written to the PDF
print(combined_plot)

dev.off()



