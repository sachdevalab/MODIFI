#!/usr/bin/env Rscript
# Combined K. pneumoniae ECE figure (new network, 14 clusters). Run with r_env Rscript.
#  a  cluster size + predicted mobility
#  b  known plasmid / host-genus reach / AMR / metal-biocide cargo (per-column-scaled heatmap)
#  c  Enterobacteriaceae relative abundance across the infant AND asthma cohorts (mean genome depth)
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork)})

KP  <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno/kp"
ENT <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs/entero_abundance_sourcedata.csv"

## ---------- panels a/b: 14 clusters ----------
clu <- read.csv(file.path(KP, "kp_clusters_fig.csv"), stringsAsFactors = FALSE) %>%
  arrange(desc(type == "virus"), n_members)   # viruses at the bottom, plasmids by size above
clu$cluster <- factor(clu$cluster, levels = clu$cluster)
clu$mobility <- factor(clu$mobility, levels = c("conjugative", "mobilizable", "non-mobilizable"))
MOB_COL <- c("conjugative" = "#cb181d", "mobilizable" = "#fd8d3c", "non-mobilizable" = "#bdbdbd")
virus_clu <- as.character(clu$cluster[clu$type == "virus"])   # label viruses on the y-axis

pa <- ggplot(clu, aes(n_members, cluster, fill = mobility)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = n_members), hjust = -0.3, size = 2.9, color = "grey20") +
  scale_fill_manual(values = MOB_COL, name = "Predicted mobility") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.13)), breaks = c(0, 10, 20, 30)) +
  scale_y_discrete(labels = function(b) ifelse(b %in% virus_clu, paste0(b, " (virus)"), paste0(b, " (plasmid)"))) +
  labs(x = "ECEs in cluster (95% ANI)", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 9), legend.position = "bottom")

tile <- clu %>% transmute(cluster,
                          `Known\nplasmid` = known,
                          `Host\ngenera`   = host_genera,
                          `AMR\ngenes`     = n_amr,
                          `Metal/\nbiocide` = n_metal) %>%
  pivot_longer(-cluster, names_to = "feature", values_to = "v")
tile$feature <- factor(tile$feature, levels = c("Known\nplasmid","Host\ngenera","AMR\ngenes","Metal/\nbiocide"))
tile$cluster <- factor(tile$cluster, levels = levels(clu$cluster))
tile <- tile %>% group_by(feature) %>%
  mutate(fill01 = if (max(v) > 0) v / max(v) else 0,
         is_known = feature == "Known\nplasmid",
         lab = ifelse(is_known, ifelse(v == 1, "yes", "–"), as.character(v))) %>% ungroup()
tile$fillval   <- ifelse(tile$is_known, NA, tile$fill01)
tile$fillknown <- ifelse(tile$is_known, tile$fill01, NA)

pb <- ggplot(tile, aes(feature, cluster)) +
  geom_tile(aes(fill = fillval), color = "white", linewidth = 0.7) +
  geom_tile(aes(fill = NULL, alpha = fillknown), fill = "#238b45", color = "white", linewidth = 0.7) +
  geom_text(aes(label = lab), size = 2.9, color = "grey15") +
  scale_fill_gradient(low = "#eaf1f8", high = "#08519c", na.value = "transparent", guide = "none") +
  scale_alpha_continuous(range = c(0, 0.9), na.value = 0, guide = "none") +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(size = 9, lineheight = 0.9), axis.ticks = element_blank())

## ---------- panel c: Enterobacteriaceae abundance, infant + asthma ----------
d0 <- read.csv(ENT, stringsAsFactors = FALSE)
clean_name <- function(s) gsub("_[A-Z]+ ", " ", gsub("_[A-Z]+$", "", s))
d0$disp <- clean_name(d0$species)
d0 <- d0 %>% group_by(sample, cohort, disp) %>%
  summarise(v = sum(rel_abundance_among_genomes), .groups = "drop")
keep <- (d0 %>% group_by(disp) %>% summarise(m = max(v), .groups = "drop") %>% filter(m >= 0.03))$disp
d0$species_lab <- ifelse(d0$disp %in% keep, d0$disp, "Other Enterobacteriaceae")
d0 <- d0 %>% group_by(sample, cohort, species_lab) %>% summarise(v = sum(v), .groups = "drop")

# sample order within each cohort (numeric)
ord_samples <- function(pre) {
  n <- sort(unique(as.integer(sub(paste0(pre, "_"), "", grep(paste0("^", pre, "_"), unique(d0$sample), value = TRUE)))))
  paste0(pre, "_", n)
}
lev <- c(ord_samples("infant"), ord_samples("asthma"))
d0$sample <- factor(d0$sample, levels = lev)

sp_all <- unique(d0$species_lab)
sp_order <- c(sort(setdiff(sp_all, "Other Enterobacteriaceae")), "Other Enterobacteriaceae")
kleb_sp  <- grep("Klebsiella", sp_order, value = TRUE)
other_sp <- setdiff(sp_order, c(kleb_sp, "Other Enterobacteriaceae"))
kleb_cols  <- colorRampPalette(c("#7f0000", "#d7301f", "#fc8d59"))(max(length(kleb_sp), 1))
other_cols <- colorRampPalette(c("#08519c","#4292c6","#41ab5d","#a1d99b","#807dba","#c994c7"))(length(other_sp))
names(kleb_cols) <- kleb_sp; names(other_cols) <- other_sp
pal <- c(kleb_cols, other_cols, "Other Enterobacteriaceae" = "#bdbdbd")
ital <- function(x) sapply(x, function(s)
  if (s == "Other Enterobacteriaceae") "'Other Enterobacteriaceae'" else paste0("italic('", s, "')"))
d0$species_lab <- factor(d0$species_lab, levels = rev(sp_order))

pc <- ggplot(d0, aes(sample, v, fill = species_lab)) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  scale_fill_manual(values = pal, breaks = sp_order,
                    labels = parse(text = ital(sp_order)), name = "Species", drop = FALSE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Enterobacteriaceae relative\nabundance") +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.text = element_text(size = 8), legend.key.size = unit(0.85, "lines"))

## ---------- combine (abundance = a on top; cluster = b, heatmap = c below) ----------
bottom <- pa + pb + plot_layout(widths = c(1.7, 1))
fig <- pc / bottom + plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(KP, "fig_kp_combined.pdf"), fig, width = 13, height = 11)
ggsave(file.path(KP, "fig_kp_combined.png"), fig, width = 13, height = 11, dpi = 300)
cat("wrote fig_kp_combined.pdf/.png\n")
