library(ggplot2)
library(readr)
library(gridExtra)
library(cowplot)


pdf(file="../../tmp/results/ecoli_base_combine.pdf", width=8, height=5, onefile=FALSE)

# Read the CSV file
df <- read_csv("../../tmp/results/ecoli_base2.csv")
df2 <- read_csv("../../tmp/results/ecoli_base_meta2.csv")

# Create the line plot with color representing coverage
p1 <- ggplot(df, aes(x = FDR, y = recall, color = as.factor(coverage))) +
  geom_line() +
  theme_minimal() +
  labs(x = "FP", y = "Recall", color = "Coverage") +
   ggtitle("C227")+
  xlim(0, 0.05)  

p2 <- ggplot(df2, aes(x = FDR, y = recall, color = as.factor(coverage))) +
  geom_line() +
  theme_minimal() +
  labs(x = "FP", y = "Recall", color = "Coverage")+
   ggtitle("C227+soil")+
  xlim(0, 0.05) 


# Combine plots without legends
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

# Explicitly print the combined plot to ensure it's written to the PDF
print(combined_plot)

dev.off()



