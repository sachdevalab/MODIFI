#!/usr/bin/env Rscript
# Panel d (PB24) recovery vs depth, split by modification type, aggregated POOLED-BY-LEVEL so the
# curves are monotonic. Fork of plot_motif_vs_depth_bytype.R. Original untouched.
#
# Why: the depth-binned mean-of-per-contig-ratios in the original wiggles because different 10x
# bins average different contig subsets (a contig whose motif only clears the 0.6/500 gate at high
# depth contributes 0% to a middle bin). Every per-contig count is actually monotonic in depth, so
# aggregating by SUBSAMPLE LEVEL with a pooled (micro-averaged) ratio is monotonic by construction.
#
# recovery(type, level) = 100 * sum_contigs count(contig,type,level) / sum_contigs count(contig,type,p100)
# x-axis: one shared mean depth per level (mean of mean_depth over all panel contigs at that level);
#         each tick dual-labelled with proportion (5%..100%) AND mean depth (Nx). 95% CI = bootstrap
#         over contigs.
#
# Reads : /home/shuaiw/borg/revision/motif_class/coverage_motif_summary_bytype.csv
# Writes: /home/shuaiw/MODIFI/tmp/rev_figs/motif_classification/
#           panel_d_by_modtype_pooled.pdf/.png  +  panel_d_by_modtype_pooled_sourcedata.csv

suppressMessages({library(ggplot2); library(dplyr)})
set.seed(42)

data_path <- "/home/shuaiw/borg/revision/motif_class/coverage_motif_summary_bytype.csv"
out_dir   <- "/home/shuaiw/MODIFI/tmp/rev_figs/motif_classification"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

plot_types <- c("6mA", "4mC", "unclassified")
type_cols  <- c("6mA" = "#C0392B", "4mC" = "#2E86C1", "unclassified" = "#E59400")
N_BOOT <- 1000

d <- read.csv(data_path, stringsAsFactors = FALSE)
d$p <- as.integer(sub(".*\\.p([0-9]+)$", "\\1", d$folder))

# Match the original panel-d contig set: drop contigs whose TOTAL motif count (summed over types,
# max over levels) exceeds 120.
tot <- d %>% group_by(folder, contig) %>% summarise(t = sum(motif_count), .groups = "drop") %>%
  group_by(contig) %>% summarise(mx = max(t), .groups = "drop")
excl <- tot %>% filter(mx > 120) %>% pull(contig)
if (length(excl)) message("Excluding high-motif contigs: ", paste(excl, collapse = ", "))
d <- d %>% filter(!contig %in% excl)

levels_p <- sort(unique(d$p))
p100 <- max(levels_p)

# ---- shared per-level mean depth (over all panel contigs present at that level) ----
lvl_depth <- d %>% distinct(contig, p, depth) %>%
  group_by(p) %>% summarise(mean_depth = mean(depth, na.rm = TRUE), .groups = "drop")

# ---- pooled recovery per (type, level) + bootstrap CI over contigs ----
pooled_one <- function(sub, contigs) {
  # sub: rows for one type; contigs: vector to include (for bootstrap resampling with dup names)
  num <- sapply(levels_p, function(L)
    sum(sub$motif_count[sub$p == L & sub$contig %in% contigs]))
  den <- sum(sub$motif_count[sub$p == p100 & sub$contig %in% contigs])
  if (den == 0) return(rep(NA_real_, length(levels_p)))
  100 * num / den
}

# bootstrap that respects duplicated contigs (weight by multiplicity)
pooled_weighted <- function(sub, tab) {
  # tab: named integer vector contig -> multiplicity
  num <- sapply(levels_p, function(L) {
    r <- sub[sub$p == L, ]
    sum(r$motif_count * ifelse(is.na(tab[r$contig]), 0, tab[r$contig]))
  })
  r100 <- sub[sub$p == p100, ]
  den <- sum(r100$motif_count * ifelse(is.na(tab[r100$contig]), 0, tab[r100$contig]))
  if (den == 0) return(rep(NA_real_, length(levels_p)))
  100 * num / den
}

res <- list()
for (ty in plot_types) {
  sub <- d %>% filter(mod_type == ty) %>% select(contig, p, motif_count)
  # contigs that carry this type at full depth (denominator basis)
  ctgs <- sub %>% filter(p == p100, motif_count > 0) %>% pull(contig) %>% unique()
  if (!length(ctgs)) next
  sub <- sub %>% filter(contig %in% ctgs)
  point <- pooled_one(sub, ctgs)

  # bootstrap over contigs
  boot <- matrix(NA_real_, nrow = N_BOOT, ncol = length(levels_p))
  for (b in seq_len(N_BOOT)) {
    samp <- sample(ctgs, length(ctgs), replace = TRUE)
    tab <- table(samp)
    boot[b, ] <- pooled_weighted(sub, tab)
  }
  lo <- apply(boot, 2, quantile, probs = 0.025, na.rm = TRUE)
  hi <- apply(boot, 2, quantile, probs = 0.975, na.rm = TRUE)

  res[[ty]] <- data.frame(mod_type = ty, p = levels_p,
                          recovery = point, lo = lo, hi = hi,
                          n_contigs = length(ctgs))
}
df <- bind_rows(res) %>% left_join(lvl_depth, by = "p") %>%
  mutate(mod_type = factor(mod_type, levels = plot_types))

# monotonicity check
for (ty in plot_types) {
  v <- df$recovery[df$mod_type == ty & order(df$p[df$mod_type == ty])]
  v <- df %>% filter(mod_type == ty) %>% arrange(p) %>% pull(recovery)
  message(sprintf("%-12s recovery p05->p100: %s  monotonic=%s  ends@%.0f%%",
                  ty, paste(round(v), collapse = ","),
                  all(diff(v) >= -1e-9), tail(v, 1)))
}

write.csv(df %>% select(mod_type, p, mean_depth, n_contigs, recovery, lo, hi),
          file.path(out_dir, "panel_d_by_modtype_pooled_sourcedata.csv"), row.names = FALSE)

# ---- dual x tick labels: proportion + mean depth ----
brk <- lvl_depth$mean_depth
lab <- paste0(lvl_depth$p, "%\n", round(lvl_depth$mean_depth), "x")

# small per-type horizontal dodge (in depth units), proportional to axis span
span <- diff(range(brk))
off <- c("6mA" = -0.008, "4mC" = 0, "unclassified" = 0.008) * span
df$x <- df$mean_depth + off[as.character(df$mod_type)]

lab_leg <- sapply(plot_types, function(t) {
  n <- unique(df$n_contigs[df$mod_type == t]); paste0(t, " (n=", n, ")")
})

p <- ggplot(df, aes(x = x, y = recovery, color = mod_type, fill = mod_type)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.6) +
  scale_color_manual("Modification type", values = type_cols,
                     breaks = plot_types, labels = lab_leg[plot_types]) +
  scale_fill_manual(guide = "none", values = type_cols) +
  scale_x_continuous("Subsample proportion / mean depth",
                     breaks = brk, labels = lab) +
  scale_y_continuous("Motif detection recovery rate (%)",
                     limits = c(0, 105), breaks = seq(0, 100, by = 20)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major = element_line(linewidth = 0.3, colour = "grey85"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 9, lineheight = 0.9),
        legend.position = c(0.98, 0.05), legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.3))

ggsave(file.path(out_dir, "panel_d_by_modtype_pooled.pdf"), p, width = 6.2, height = 3.9)
ggsave(file.path(out_dir, "panel_d_by_modtype_pooled.png"), p, width = 6.2, height = 3.9, dpi = 200)
message("Saved: ", file.path(out_dir, "panel_d_by_modtype_pooled.pdf"))
