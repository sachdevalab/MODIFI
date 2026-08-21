#!/usr/bin/env Rscript
# Step 4 figure - Enterobacteriaceae abundance per infant gut (reviewer point 5.1), two panels:
#  (a) composition WITHIN the Enterobacteriaceae community (bars sum to 100%) - is Klebsiella
#      merely the dominant Enterobacteriaceae?
#  (b) ABSOLUTE Enterobacteriaceae abundance in the WHOLE metagenome (bars = Enterobacteriaceae
#      fraction of the classified community) - controls for the concern that a taxon can dominate
#      Enterobacteriaceae yet be low-biomass overall (few ECEs to detect).
# Shared species colours + one legend. Species names italicised. Run with r_env Rscript.

suppressMessages({
  library(ggplot2); library(dplyr); library(scales); library(patchwork)
})

fig_dir <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
d0 <- read.csv(file.path(fig_dir, "entero_abundance_sourcedata.csv"), stringsAsFactors = FALSE)

# clean display name: drop GTDB alphabetic suffixes (Serratia_J -> Serratia, marcescens_K -> ...)
clean_name <- function(s) gsub("_[A-Z]+ ", " ", gsub("_[A-Z]+$", "", s))
d0$disp <- clean_name(d0$species); d0$genus <- clean_name(d0$genus)
d0 <- d0 %>% group_by(sample, disp, genus) %>%
  summarise(rel_abundance_entero = sum(rel_abundance_entero),
            rel_abundance_among_genomes = sum(rel_abundance_among_genomes), .groups = "drop")

# species kept in the legend: >=3% within Enterobacteriaceae in some gut; else "Other"
maxrel <- d0 %>% group_by(disp) %>% summarise(m = max(rel_abundance_entero), .groups = "drop")
keep <- maxrel$disp[maxrel$m >= 0.03]
d0$species_lab <- ifelse(d0$disp %in% keep, d0$disp, "Other Enterobacteriaceae")

# sample order (infant_1 .. infant_28)
snum <- as.integer(sub("infant_", "", d0$sample))
lev <- paste0("infant_", sort(unique(snum)))

agg <- function(df, val) df %>% group_by(sample, species_lab, genus) %>%
  summarise(v = sum(.data[[val]]), .groups = "drop") %>%
  mutate(sample = factor(sample, levels = lev))

# consistent species ordering + palette (Klebsiella reds, others cool colours)
sp_all <- d0 %>% distinct(species_lab) %>% pull(species_lab)
sp_order <- c(sort(setdiff(sp_all, "Other Enterobacteriaceae")), "Other Enterobacteriaceae")
kleb_sp  <- grep("Klebsiella", sp_order, value = TRUE)
other_sp <- setdiff(sp_order, c(kleb_sp, "Other Enterobacteriaceae"))
kleb_cols  <- colorRampPalette(c("#7f0000", "#d7301f", "#fc8d59"))(max(length(kleb_sp), 1))
other_cols <- colorRampPalette(c("#08519c", "#4292c6", "#41ab5d", "#a1d99b",
                                 "#807dba", "#c994c7"))(length(other_sp))
names(kleb_cols) <- kleb_sp; names(other_cols) <- other_sp
pal <- c(kleb_cols, other_cols, "Other Enterobacteriaceae" = "#bdbdbd")

ital <- function(x) sapply(x, function(s)
  if (s == "Other Enterobacteriaceae") "'Other Enterobacteriaceae'" else paste0("italic('", s, "')"))

mk <- function(df, ytitle, title, pct = TRUE) {
  df$species_lab <- factor(df$species_lab, levels = rev(sp_order))
  ggplot(df, aes(sample, v, fill = species_lab)) +
    geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
    scale_fill_manual(values = pal, breaks = sp_order,
                      labels = parse(text = ital(sp_order)), name = "Species", drop = FALSE) +
    scale_y_continuous(labels = if (pct) percent_format(accuracy = 1) else waiver(),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = ytitle, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(size = 8), legend.key.size = unit(0.9, "lines"),
          plot.title = element_text(size = 11, face = "plain"))
}

pb <- mk(agg(d0, "rel_abundance_among_genomes"),
         "Enterobacteriaceae relative abundance (mean depth)",
         "Enterobacteriaceae abundance in infant guts (mean genome depth, fraction of all reconstructed genomes)")

ggsave(file.path(fig_dir, "fig_entero_abundance.pdf"), pb, width = 11, height = 5.2)
ggsave(file.path(fig_dir, "fig_entero_abundance.png"), pb, width = 11, height = 5.2, dpi = 300)
cat("wrote fig_entero_abundance.pdf/.png (panel b only)\n")
