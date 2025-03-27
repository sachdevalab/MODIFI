
# Load necessary libraries
library(ggplot2)

# Load the CSV file into a DataFrame
file_path <- "/home/shuaiw/methylation/data/borg/human/test_result6.csv"
df <- read.csv(file_path)

# Print the first few rows of the DataFrame to verify
print(head(df))


## calculate the correlation coefficient and makr the plot
# Calculate the Pearson correlation coefficient
correlation <- cor(df$estimated, df$real, method = "pearson")
correlation_text <- paste("Pearson correlation: ", round(correlation, 2))


# Create a scatter plot using ggplot2
p1 <- ggplot(df, aes(x = real, y = estimated)) +
  geom_point(size = 2) +
  labs(x = "Observed", y = "Predicted") +
  # ggtitle("Scatter Plot of Predicted vs Observed IPD") +

  # theme_minimal()+
  annotate("text", x = Inf, y = Inf, label = correlation_text, hjust = 1.1, vjust = 2, size = 8, color = "blue")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlim(0, 5) +
  ylim(0, 5)


# Calculate the Pearson correlation coefficient
correlation <- cor(df$ipd_sum, df$real, method = "pearson")
correlation_text <- paste("Pearson correlation: ", round(correlation, 2))


# Create a scatter plot using ggplot2
p2 <- ggplot(df, aes(x = real, y = ipd_sum)) +
  geom_point(size = 2) +
  labs(x = "Observed", y = "Predicted") +
  # ggtitle("Scatter Plot of Predicted vs Observed IPD") +

  # theme_minimal()+
  annotate("text", x = Inf, y = Inf, label = correlation_text, hjust = 1.1, vjust = 2, size = 8, color = "blue")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  xlim(0, 5) +
  ylim(0, 5)

## plot p1 and p2 together
library(gridExtra)
combined_plot <- grid.arrange(p1, p2, ncol=2)




# Save the plot to a file
ggsave("tmp/scatter_plot_meta.pdf", combined_plot, width = 12, height = 6)