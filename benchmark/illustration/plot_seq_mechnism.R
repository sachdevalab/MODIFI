# Required packages
library(ggplot2)
library(dplyr)

# Define simulation parameters
set.seed(42)
bases <- c("T", "C", "G", "A", "T", "C", "A", "A")
base_colors <- c("T" = "green", "C" = "orange", "G" = "red", "A" = "blue")
base_times_mod <- c(70.0, 70.4, 71.2, 71.9, 73.5, 73.8, 74.0, 74.3)
base_ctrl <- c("T", "C", "G", "A", "T", "C", "A", "A", "G", "C", "T", "A")
base_times_ctrl <- c(70.0, 70.4, 71.2, 71.9, 72.4, 72.7, 72.9, 73.2, 73.5, 73.9, 74.2, 74.5)
duration <- 0.3
sampling_rate <- 0.01
time <- seq(69.8, 74.8, by = sampling_rate)

pdf("./plot_seq_mechnism.pdf", width = 8, height = 5)

simulate_trace <- function(base_times, noise_sd = 5) {
  signal <- rnorm(length(time), mean = 10, sd = noise_sd)
  for (bt in base_times) {
    idx <- which.min(abs(time - bt))
    width <- round(duration / sampling_rate)
    peak <- dnorm(seq(-2, 2, length.out = width)) * 300
    start <- max(1, idx - floor(width/2))
    end <- min(length(time), start + width - 1)
    signal[start:end] <- signal[start:end] + peak[1:(end - start + 1)]
  }
  return(signal)
}

# Simulate both traces
signal_mod <- simulate_trace(base_times_mod)
signal_ctrl <- simulate_trace(base_times_ctrl)

# Create data frame for ggplot
df_mod <- data.frame(time = time, signal = signal_mod)
df_ctrl <- data.frame(time = time, signal = signal_ctrl)

# Pulse labels for modification
labels <- data.frame(
  base = bases,
  time_mod = base_times_mod,
  color = base_colors[bases]
)

# Pulse labels for control
labels_ctrl <- data.frame(
  base_ctrl = base_ctrl,
  time_ctrl = base_times_ctrl,
  color = base_colors[base_ctrl]
)

# Plot with modification
p1<-ggplot(df_mod, aes(x = time, y = signal)) +
  geom_line(color = "gray50", alpha = 0.6) +
  lapply(1:nrow(labels), function(i) {
    geom_line(
      data = df_mod %>% filter(time >= labels$time_mod[i] - 0.15 & time <= labels$time_mod[i] + 0.15),
      aes(x = time, y = signal), color = labels$color[i], linewidth = 1
    )
  }) +
  geom_text(data = labels, aes(x = time_mod, y = 220, label = base, color = base), size = 5, show.legend = FALSE) +
  scale_color_manual(values = base_colors) +
  annotate("segment", x = 71.9, xend = 73.5, y = 360, yend = 360, linetype = "dashed", arrow = arrow(length = unit(0.2,"cm"))) +
  annotate("text", x = 72.7, y = 365, label = "mA") +
  labs(x = "Time (s)", y = "Fluorescence Intensity (a.u.)") +
  ylim(0, 250) +
  theme_minimal()

# Plot without modification (control)
p2<-ggplot(df_ctrl, aes(x = time, y = signal)) +
  geom_line(color = "gray50", alpha = 0.6) +
  lapply(1:nrow(labels_ctrl), function(i) {
    geom_line(
      data = df_ctrl %>% filter(time >= labels_ctrl$time_ctrl[i] - 0.15 & time <= labels_ctrl$time_ctrl[i] + 0.15),
      aes(x = time, y = signal), color = labels_ctrl$color[i], linewidth = 1
    )
  }) +
  geom_text(data = labels_ctrl, aes(x = time_ctrl, y = 220, label = base_ctrl, color = base_ctrl), size = 5, show.legend = FALSE) +
  scale_color_manual(values = base_colors) +
  labs(x = "Time (s)", y = "Fluorescence Intensity (a.u.)") +
  ylim(0, 250) +
  theme_minimal()

# Combine plots
library(gridExtra)
combined_plot <- gridExtra::grid.arrange(p1, p2, ncol = 1)

## save to pdf
# ggsave("./plot_seq_mechnism.pdf", width = 10, height = 5)

