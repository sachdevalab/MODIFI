library(ggplot2)
library(dplyr)

# Hi-C validation boxplot for the STRICT ECE dataset (R version of plot_bin3c.R).
# Reads the strict data in tmp/rev_figs/ece_anno/hic/ (category column "ECE type", already labeled
# "MODIFI linkages supported by Hi-C" vs "Random Pairs"). Legend title = "Method".
fig_dir <- "../../tmp/rev_figs/ece_anno/hic/"
df <- read.csv(file.path(fig_dir, "contact_values.csv"))
LAB <- "MODIFI linkages supported by Hi-C"

our_v <- df %>% filter(ECE.type == LAB) %>% pull(Contact.Value)
rand_v <- df %>% filter(ECE.type == "Random Pairs") %>% pull(Contact.Value)
tt <- t.test(our_v, rand_v, var.equal = FALSE)
mw <- wilcox.test(our_v, rand_v, alternative = "greater")
cat(sprintf("n linkages=%d  mean=%.1f median=%.0f | random mean=%.3f\n",
            length(our_v), mean(our_v), median(our_v), mean(rand_v)))
cat(sprintf("Welch t p=%.3e ; Mann-Whitney U p=%.3e ; positive %d/%d\n",
            tt$p.value, mw$p.value, sum(our_v > 0), length(our_v)))

# +1 pseudocount so zero-contact random pairs still render on the log axis
df$Contact.Value.Log <- df$Contact.Value + 1
df$ECE.type <- factor(df$ECE.type, levels = c(LAB, "Random Pairs"))

p <- ggplot(df, aes(x = Sample, y = Contact.Value.Log, fill = ECE.type)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(y = "Contact Value + 1 (log scale)", x = "", fill = "Method", title = "") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top", legend.text = element_text(size = 8),
        legend.title = element_text(size = 9))

ggsave(file.path(fig_dir, "bin3c_boxplot.pdf"), p, width = 4.5, height = 5.3)
ggsave(file.path(fig_dir, "bin3c_boxplot.png"), p, width = 4.5, height = 5.3, dpi = 200)
cat("wrote", file.path(fig_dir, "bin3c_boxplot.{pdf,png}"), "\n")
