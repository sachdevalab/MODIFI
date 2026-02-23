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

# Find maximum values and print them
max_cpu_idx <- which.max(df$cpu_time_hours)
max_wallclock_idx <- which.max(df$wallclock_time_hours)
max_memory_idx <- which.max(df$peak_memory_gb)

cat("=== CPU Time Statistics ===\n")
cat("Min: ", min(df$cpu_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Max: ", max(df$cpu_time_hours, na.rm = TRUE), " hours (Sample: ", df$sample[max_cpu_idx], ", Env: ", df$environment[max_cpu_idx], ")\n", sep="")
cat("Mean: ", mean(df$cpu_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Std: ", sd(df$cpu_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Range: ", max(df$cpu_time_hours, na.rm = TRUE) - min(df$cpu_time_hours, na.rm = TRUE), " hours\n\n", sep="")

cat("=== Wall-clock Time Statistics ===\n")
cat("Min: ", min(df$wallclock_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Max: ", max(df$wallclock_time_hours, na.rm = TRUE), " hours (Sample: ", df$sample[max_wallclock_idx], ", Env: ", df$environment[max_wallclock_idx], ")\n", sep="")
cat("Mean: ", mean(df$wallclock_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Std: ", sd(df$wallclock_time_hours, na.rm = TRUE), " hours\n", sep="")
cat("Range: ", max(df$wallclock_time_hours, na.rm = TRUE) - min(df$wallclock_time_hours, na.rm = TRUE), " hours\n\n", sep="")

cat("=== Peak Memory Statistics ===\n")
cat("Min: ", min(df$peak_memory_gb, na.rm = TRUE), " GB\n", sep="")
cat("Max: ", max(df$peak_memory_gb, na.rm = TRUE), " GB (Sample: ", df$sample[max_memory_idx], ", Env: ", df$environment[max_memory_idx], ")\n", sep="")
cat("Mean: ", mean(df$peak_memory_gb, na.rm = TRUE), " GB\n", sep="")
cat("Std: ", sd(df$peak_memory_gb, na.rm = TRUE), " GB\n", sep="")
cat("Range: ", max(df$peak_memory_gb, na.rm = TRUE) - min(df$peak_memory_gb, na.rm = TRUE), " GB\n\n", sep="")

# Create CPU time box plot
max_cpu_data <- df[max_cpu_idx, ]
p1 <- ggplot(df, aes(x = environment, y = cpu_time_hours, fill = environment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  # geom_point(data = max_cpu_data, aes(x = environment, y = cpu_time_hours), 
  #            color = "red", size = 4, shape = 18) +
  # geom_text(data = max_cpu_data, aes(x = environment, y = cpu_time_hours, 
  #                                     label = sample),
  #           vjust = -1, hjust = 0.5, color = "red", size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       x = "Environment",
       y = "CPU Time (hours, log scale)") +
  scale_y_log10() +
  scale_fill_manual(values = env_colors) +
  coord_cartesian(clip = "off", ylim = c(NA, max(df$cpu_time_hours) * 2.5))

# Create wall-clock time box plot
max_wallclock_data <- df[max_wallclock_idx, ]
p2 <- ggplot(df, aes(x = environment, y = wallclock_time_hours, fill = environment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  # geom_point(data = max_wallclock_data, aes(x = environment, y = wallclock_time_hours), 
  #            color = "red", size = 4, shape = 18) +
  # geom_text(data = max_wallclock_data, aes(x = environment, y = wallclock_time_hours, 
  #                                           label = sample),
  #           vjust = -1, hjust = 0.5, color = "red", size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       x = "Environment",
       y = "Wall-clock Time (hours, log scale)") +
  scale_y_log10() +
  scale_fill_manual(values = env_colors) +
  coord_cartesian(clip = "off", ylim = c(NA, max(df$wallclock_time_hours) * 2.5))

# Create peak memory box plot
max_memory_data <- df[max_memory_idx, ]
p3 <- ggplot(df, aes(x = environment, y = peak_memory_gb, fill = environment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  # geom_point(data = max_memory_data, aes(x = environment, y = peak_memory_gb), 
  #            color = "red", size = 4, shape = 18) +
  # geom_text(data = max_memory_data, aes(x = environment, y = peak_memory_gb, 
  #                                        label = sample),
  #           vjust = -1, hjust = 0.5, color = "red", size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(title = "",
       x = "Environment",
       y = "Peak Memory (GB, log scale)") +
  scale_y_log10() +
  scale_fill_manual(values = env_colors) +
  coord_cartesian(clip = "off", ylim = c(NA, max(df$peak_memory_gb) * 2.5))

# Combine plots
combined_plot <- p1 + p2 + p3

# Display the plot
print(combined_plot)

# Save the plot
ggsave("/home/shuaiw/mGlu/tmp/figures/multi_env_linkage/time_usage_by_environment.pdf", 
       combined_plot, 
       width = 12, 
       height = 4)
