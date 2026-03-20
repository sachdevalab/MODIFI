# ipdSummary vs MODIFI: CPU, wall time, and peak RAM vs reference size
#
# Run from anywhere:
#   conda run -n r_env Rscript /home/shuaiw/MODIFI/benchmark/compare/plot_runtime_benchmark.R
#
# Requires: tidyverse (readr, ggplot2, dplyr, tidyr, scales)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

modifi_root <- Sys.getenv("MODIFI_ROOT", unset = "/home/shuaiw/MODIFI")
in_csv <- file.path(modifi_root, "tmp/figures/base_benchmark/ipd_vs_modifi_runtime_long.csv")
out_pdf <- file.path(modifi_root, "tmp/figures/base_benchmark/runtime_vs_refsize.pdf")
out_png <- file.path(modifi_root, "tmp/figures/base_benchmark/runtime_vs_refsize.png")

stopifnot(file.exists(in_csv))

raw <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(
    ref_mb = ref_bases / 1e6,
    software = factor(software, levels = c("MODIFI", "ipdSummary"))
  )

df <- raw %>%
  pivot_longer(
    cols = c(cpu_hr, wall_hr, max_rss_gb),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("cpu_hr", "wall_hr", "max_rss_gb"),
      labels = c(
        "CPU time (h)",
        "Wall-clock time (h)",
        "Peak Memory (GB)"
      )
    )
  ) %>%
  filter(!is.na(value))

pal <- c(MODIFI = "#0f766e", ipdSummary = "#b45309")

p <- ggplot(df, aes(x = ref_mb, y = value, color = software, group = software)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) +
  scale_color_manual(values = pal, name = NULL) +
  scale_x_continuous(
    name = "Reference size (Mb)",
    breaks = pretty_breaks(n = 6),
    labels = label_number(accuracy = 1)
  ) +
  scale_y_continuous(
    name = NULL,
    labels = label_number(accuracy = 0.1)
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(out_pdf, p, width = 11, height = 4.2, dpi = 300)
ggsave(out_png, p, width = 11, height = 4.2, dpi = 300)

message("Wrote ", normalizePath(out_pdf, winslash = "/", mustWork = FALSE))
message("Wrote ", normalizePath(out_png, winslash = "/", mustWork = FALSE))
