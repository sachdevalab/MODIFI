library(ggplot2)
library(dplyr)

# Set paths
fig_dir <- "../../tmp/figures/link_accuracy/"
df <- read.csv(file.path(fig_dir, "contact_values.csv"))

# T-test results
our_linkages_values <- df %>% filter(Type == "Our Linkages") %>% pull(Contact.Value)
random_pairs_values <- df %>% filter(Type == "Random Pairs") %>% pull(Contact.Value)
t_test_result <- t.test(our_linkages_values, random_pairs_values, var.equal = FALSE)
cat(sprintf("T-test results: t-statistic = %.4f, p-value = %.4e\n", 
            t_test_result$statistic, t_test_result$p.value))
cat(sprintf("Average contact value for our linkages: %.4f\n", mean(our_linkages_values)))
cat(sprintf("Average contact value for random pairs: %.4f\n", mean(random_pairs_values)))

# T-test for each sample
for (sample in unique(df$Sample)) {
    sample_our_linkages <- df %>% filter(Type == "Our Linkages", Sample == sample) %>% pull(Contact.Value)
    sample_random_pairs <- df %>% filter(Type == "Random Pairs", Sample == sample) %>% pull(Contact.Value)
    sample_t_test <- t.test(sample_our_linkages, sample_random_pairs, var.equal = FALSE)
    cat(sprintf("Sample %s T-test results: t-statistic = %.4f, p-value = %.4e\n",
                sample, sample_t_test$statistic, sample_t_test$p.value))
}

# Print dataframe summary
# print(df)
cat(sprintf("\nAverage contact value for our linkages: %.4f\n", 
            mean(df %>% filter(Type == "Our Linkages") %>% pull(Contact.Value))))
cat(sprintf("Average contact value for random pairs: %.4f\n", 
            mean(df %>% filter(Type == "Random Pairs") %>% pull(Contact.Value))))

# Create plot with pseudocount to handle zeros
df$Contact.Value.Log <- df$Contact.Value + 1  # Add 1 before log transformation

p <- ggplot(df, aes(x = Sample, y = Contact.Value.Log, fill = Type)) +
    geom_boxplot() +
    scale_y_log10() +
    labs(y = "Contact Value + 1 (log scale)", title = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "top")

ggsave(file.path(fig_dir, "bin3c_boxplot.pdf"), p, width = 4, height = 5)
