library(ggplot2)
library(dplyr)
library(patchwork)

fig_dir   <- "../../tmp/figures/link_accuracy/"
df        <- read.csv(file.path(fig_dir, "linkage_curves_all.csv"))
type_df   <- read.csv(file.path(fig_dir, "ece_type_counts.csv"))

# Drop the origin row (score_threshold = NA = no cutoff applied)
df <- df %>% filter(!is.na(score_threshold))

df <- df %>%
  mutate(n_species = factor(n_species, levels = c(10, 20, 30, 40, 50)))

# ── Palette ────────────────────────────────────────────────────────────────────
palette <- c(
  "10" = "#1B7837",
  "20" = "#4393C3",
  "30" = "#D6604D",
  "40" = "#8073AC",
  "50" = "#E08214"
)

# ── Theme ──────────────────────────────────────────────────────────────────────
base_theme <- theme_classic(base_size = 12) +
  theme(
    panel.grid.major = element_line(colour = "grey93", linewidth = 0.35),
    axis.line        = element_line(colour = "grey40"),
    axis.ticks       = element_line(colour = "grey40"),
    axis.title       = element_text(size = 11),
    legend.title     = element_text(size = 10, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.width = unit(1.5, "lines"),
    legend.position  = "right"
  )

# ── Panel A: Recall vs score cutoff ───────────────────────────────────────────
p_recall <- ggplot(df, aes(x = score_threshold, y = recall,
                            colour = n_species, group = sample)) +
  geom_step(direction = "vh", linewidth = 0.9) +
  scale_colour_manual(values = palette, name = "Community\nsize") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(add = 0.01)) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(add = 0)
  ) +
  labs(x = "Score cutoff", y = "Recall") +
  base_theme

# ── Panel B: Precision vs score cutoff ────────────────────────────────────────
p_precision <- ggplot(df, aes(x = score_threshold, y = precision,
                               colour = n_species, group = sample)) +
  geom_step(direction = "vh", linewidth = 0.9) +
  scale_colour_manual(values = palette, name = "Community\nsize") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(add = 0.01)) +
  scale_y_continuous(
    limits = c(0, 1.05),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(add = 0)
  ) +
  labs(x = "Score cutoff", y = "Precision") +
  base_theme

# ── Panel C: ECE type barplot ──────────────────────────────────────────────────
type_df <- type_df %>%
  mutate(
    n_species = factor(n_species, levels = c(10, 20, 30, 40, 50)),
    mge_type  = factor(mge_type, levels = c("plasmid", "virus"))
  )

type_palette <- c("plasmid" = "#5B8DB8", "virus" = "#E07B54")

p_bar <- ggplot(type_df, aes(x = n_species, y = count, fill = mge_type)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = type_palette, name = "ECE type") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Community size (species)", y = "Number of ECEs") +
  base_theme +
  theme(legend.position = "right")

# ── Save ECE barplot independently ────────────────────────────────────────────
ggsave(file.path(fig_dir, "ece_type_counts.pdf"),
       p_bar, width = 6, height = 4.2)

# ── Save recall + precision combined ──────────────────────────────────────────
p_combined <- (p_recall | p_precision) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(file.path(fig_dir, "linkage_multi_score.pdf"),
       p_combined, width = 9, height = 4.2)

cat("Saved to", fig_dir, "\n")
