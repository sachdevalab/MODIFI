library(tidyverse)
library(patchwork)

# Read data
df <- read_csv("../../tmp/figures/inversion/inversion_ratio_results.csv")

# Prepare data for pie charts - reshape to long format
df_long <- df %>%
  select(sample, week, raw_crossing, inv_crossing) %>%
  pivot_longer(cols = c(raw_crossing, inv_crossing),
               names_to = "type",
               values_to = "count") %>%
  mutate(type = ifelse(type == "raw_crossing", "G1", "G2"),
         type = factor(type, levels = c("G1", "G2")))

# Calculate percentages for each sample
df_long <- df_long %>%
  group_by(sample, week) %>%
  mutate(percentage = count / sum(count) * 100,
         label = paste0(type, "\n", count, " (", round(percentage, 1), "%)"))

# Create individual pie charts for each time point
plots <- list()

for (i in 1:nrow(df)) {
  sample_data <- df_long %>% filter(sample == df$sample[i])
  
  p <- ggplot(sample_data, aes(x = "", y = count, fill = type)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1) +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("G1" = "#FF1493", "G2" = "#00BCD4")) +
    labs(fill = "Genotype", title = df$sample[i]) +
    theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12, vjust = 2)) +
    geom_text(aes(label = paste0(round(percentage, 1), "%")),
              position = position_stack(vjust = 0.5),
              size = 5, fontface = "bold", color = "white")
  
  plots[[i]] <- p
}

# Arrange plots in a grid with a single shared legend
combined_plot <- wrap_plots(plots, ncol = 4, guides = "collect") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Save plot
ggsave("../../tmp/figures/inversion/inversion_ratio_pie_charts.pdf", 
       combined_plot, width = 8, height = 2, dpi = 300)
ggsave("../../tmp/figures/inversion/inversion_ratio_pie_charts.png", 
       combined_plot, width = 14, height = 4, dpi = 300)

print(combined_plot)