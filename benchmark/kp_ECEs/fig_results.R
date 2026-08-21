#!/usr/bin/env Rscript
# Result figures for the Klebsiella-ECE reviewer response. Run with r_env Rscript.
#   Fig 1: overview of the 22 clusters (size, mobility, known-vs-novel, AMR/metal cargo)
#   Fig 2: cluster x host-genus reach (cross-Enterobacteriaceae transfer)
#   Fig 3: modification-motif sharing between Klebsiella and co-resident Enterobacteriaceae
#   Fig 4: infant_15_35_C deep dive (cross-genus members + gene cargo)
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork)})

D <- "/home/shuaiw/MODIFI/tmp/rev_figs/kp_ECEs"
ital <- function(x) parse(text = ifelse(grepl("^[A-Z][a-z]+ [a-z]", x),
                                        paste0("italic('", x, "')"), paste0("'", x, "'")))

GENUS_ORDER <- c("Klebsiella","Citrobacter","Enterobacter","Escherichia",
                 "Serratia","Pluralibacter","Leclercia","Salmonella")
GENUS_COL <- setNames(c("#d7301f","#2171b5","#238b45","#88419d",
                        "#f16913","#7a0177","#807dba","#525252"), GENUS_ORDER)

## ---------- load ----------
clu <- read.csv(file.path(D,"kp_22_clusters_annotated.tsv"), sep="\t", stringsAsFactors=FALSE)
clu <- clu %>% arrange(n_cluster_members)
clu$cluster_rep <- factor(clu$cluster_rep, levels=clu$cluster_rep)
# dominant mobility: conjugative > mobilizable > non-mobilizable. Note "(^|;)mobilizable" so the
# substring inside "non-mobilizable" does not wrongly match.
clu$mobility <- ifelse(grepl("conjugative", clu$predicted_mobility), "conjugative",
                ifelse(grepl("(^|;)mobilizable", clu$predicted_mobility), "mobilizable",
                       "non-mobilizable"))
clu$known <- ifelse(!is.na(suppressWarnings(as.numeric(clu$rep_mash_distance))) &
                    as.numeric(clu$rep_mash_distance) < 0.05, "known (<0.05)", "divergent / novel")

MOB_COL <- c("conjugative"="#cb181d","mobilizable"="#fd8d3c","non-mobilizable"="#bdbdbd")

## ---------- Fig 1: cluster overview ----------
clu$mobility <- factor(clu$mobility, levels=c("conjugative","mobilizable","non-mobilizable"))
p1a <- ggplot(clu, aes(n_cluster_members, cluster_rep, fill=mobility)) +
  geom_col(width=0.72) +
  geom_text(aes(label=n_cluster_members), hjust=-0.3, size=2.9, color="grey20") +
  scale_fill_manual(values=MOB_COL, name="Predicted\nmobility") +
  scale_x_continuous(expand=expansion(mult=c(0,0.13)), breaks=c(0,10,20,30)) +
  labs(x="ECEs in cluster (95% ANI)", y=NULL, title="a  Cluster size and mobility") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=11, face="bold"),
        axis.text.y=element_text(size=9), legend.position="bottom")

# panel b: heatmap with each column normalised to its own scale (so a 4-genus cell and a 43-gene
# cell are both readable); the categorical "Known plasmid" column is drawn in a separate hue.
tile <- clu %>% transmute(cluster_rep,
                          `Known\nplasmid`=ifelse(known=="known (<0.05)",1,0),
                          `Host\ngenera`=n_distinct_host_genera,
                          `AMR\ngenes`=n_amr_genes_cluster,
                          `Metal/\nbiocide`=n_stress_genes_cluster) %>%
  pivot_longer(-cluster_rep, names_to="feature", values_to="v")
tile$feature <- factor(tile$feature, levels=c("Known\nplasmid","Host\ngenera","AMR\ngenes","Metal/\nbiocide"))
tile$cluster_rep <- factor(tile$cluster_rep, levels=levels(clu$cluster_rep))
tile <- tile %>% group_by(feature) %>%
  mutate(fill01 = if (max(v)>0) v/max(v) else 0,
         is_known = feature=="Known\nplasmid",
         lab = ifelse(is_known, ifelse(v==1,"yes","–"), as.character(v))) %>% ungroup()
# two colour ramps: green for the categorical known column, blue for the count columns
tile$fillval <- ifelse(tile$is_known, NA, tile$fill01)
tile$fillknown <- ifelse(tile$is_known, tile$fill01, NA)
p1b <- ggplot(tile, aes(feature, cluster_rep)) +
  geom_tile(aes(fill=fillval), color="white", linewidth=0.7) +
  geom_tile(aes(fill=NULL, alpha=fillknown), fill="#238b45", color="white", linewidth=0.7) +
  geom_text(aes(label=lab), size=2.9, color="grey15") +
  scale_fill_gradient(low="#eaf1f8", high="#08519c", na.value="transparent", guide="none") +
  scale_alpha_continuous(range=c(0,0.9), na.value=0, guide="none") +
  scale_x_discrete(position="top") +
  labs(x=NULL,y=NULL,title="b  Known plasmid, host reach, resistance") +
  theme_minimal(base_size=11) +
  theme(axis.text.y=element_blank(), panel.grid=element_blank(),
        plot.title=element_text(size=11, face="bold"),
        axis.text.x=element_text(size=9, lineheight=0.9),
        axis.ticks=element_blank())
f1 <- p1a + p1b + plot_layout(widths=c(1.7,1))
ggsave(file.path(D,"fig_kp_clusters_overview.pdf"), f1, width=11, height=6.2)
ggsave(file.path(D,"fig_kp_clusters_overview.png"), f1, width=11, height=6.2, dpi=300)

## ---------- Fig 2: cluster x host-genus reach ----------
cg <- read.csv(file.path(D,"fig_cluster_genus_long.csv"), stringsAsFactors=FALSE)
ord <- clu %>% arrange(n_distinct_host_genera, n_cluster_members) %>% pull(cluster_rep) %>% as.character()
cg$cluster_rep <- factor(cg$cluster_rep, levels=ord)
cg$genus <- factor(cg$genus, levels=GENUS_ORDER)
cg$present <- 1
f2 <- ggplot(cg, aes(genus, cluster_rep)) +
  geom_point(aes(color=genus, size=n_members)) +
  scale_color_manual(values=GENUS_COL, guide="none") +
  scale_size_continuous(range=c(2,7), name="ECEs in cluster") +
  labs(x=NULL, y=NULL,
       title="Host-genus reach of the 22 K. pneumoniae-associated ECE clusters",
       subtitle="6 clusters (top) link across >1 Enterobacteriaceae genus; 16 are Klebsiella-restricted") +
  theme_bw(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1,face="italic"),
        plot.title=element_text(size=11), plot.subtitle=element_text(size=9),
        panel.grid.minor=element_blank())
ggsave(file.path(D,"fig_kp_crossgenus_reach.pdf"), f2, width=7.5, height=6)
ggsave(file.path(D,"fig_kp_crossgenus_reach.png"), f2, width=7.5, height=6, dpi=300)

## ---------- Fig 3: motif sharing between Klebsiella and co-resident Enterobacteriaceae ----------
ms <- read.csv(file.path(D,"kp_host_motif_similarity.tsv"), sep="\t", stringsAsFactors=FALSE)
ms$lab <- paste0(ms$sample, ": ", ms$other_species, " (", ms$other_host, ")")
ms <- ms %>% arrange(cosine)
ms$lab <- factor(ms$lab, levels=unique(ms$lab))
ms$other_genus <- factor(ms$other_genus, levels=GENUS_ORDER)
f3 <- ggplot(ms, aes(cosine, lab, color=other_genus)) +
  geom_segment(aes(x=0, xend=cosine, yend=lab), color="grey70") +
  geom_point(size=4) +
  scale_color_manual(values=GENUS_COL, name="Co-resident genus", drop=TRUE) +
  scale_x_continuous(limits=c(0,1), expand=expansion(mult=c(0,0.02))) +
  labs(x="Modification-motif similarity to co-resident Klebsiella (cosine)", y=NULL,
       title="Do Klebsiella plasmids share motifs with co-resident Enterobacteriaceae?",
       subtitle="High where cross-genus transfer occurs; low where host-restricted") +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(size=11), plot.subtitle=element_text(size=9),
        axis.text.y=element_text(size=8))
ggsave(file.path(D,"fig_kp_motif_sharing.pdf"), f3, width=9.5, height=4.6)
ggsave(file.path(D,"fig_kp_motif_sharing.png"), f3, width=9.5, height=4.6, dpi=300)

## ---------- Fig 4: infant_15_35_C deep dive ----------
mem <- read.csv(file.path(D,"infant_15_35_C_members.tsv"), sep="\t", stringsAsFactors=FALSE)
mem$genus <- ifelse(mem$network_linked_host_genus=="" | is.na(mem$network_linked_host_genus),
                    "not linked", mem$network_linked_host_genus)
gc <- c(GENUS_COL,"not linked"="#cccccc")
present <- intersect(c(GENUS_ORDER,"not linked"), unique(mem$genus))  # drop unused legend entries
mem$genus <- factor(mem$genus, levels=present)
gc <- gc[present]
mem$pid <- suppressWarnings(as.numeric(mem$pid_to_representative))
mem$pid[is.na(mem$pid)] <- 100          # the representative itself
p4a <- ggplot(mem, aes(length_bp/1000, pid, color=genus)) +
  geom_point(aes(size=aln_cov_to_representative), alpha=0.9) +
  scale_color_manual(values=gc, name="Linked host genus", drop=TRUE) +
  scale_size_continuous(range=c(2,6), name="aln cov to rep (%)") +
  scale_x_log10(limits=c(30,500)) +
  labs(x="ECE length (kb, log scale)", y="% identity to 416 kb representative",
       title="a  infant_15_35_C cluster: 24 members across 4 genera") +
  theme_classic(base_size=10) + theme(plot.title=element_text(size=10))

gc4 <- read.csv(file.path(D,"infant_15_35_C_gene_categories.csv"), stringsAsFactors=FALSE)
gc4 <- gc4 %>% filter(category!="Other / hypothetical", n_genes>0) %>% arrange(n_genes)
gc4$category <- factor(gc4$category, levels=gc4$category)
p4b <- ggplot(gc4, aes(n_genes, category)) +
  geom_col(fill="#4292c6", width=0.7) +
  geom_text(aes(label=n_genes), hjust=-0.3, size=2.8) +
  scale_x_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="genes on the 416 kb representative", y=NULL,
       title="b  infant_15_35_C functional cargo") +
  theme_classic(base_size=10) + theme(plot.title=element_text(size=10))
f4 <- p4a / p4b + plot_layout(heights=c(1.1,1))
ggsave(file.path(D,"fig_infant_15_35_C.pdf"), f4, width=8.5, height=8)
ggsave(file.path(D,"fig_infant_15_35_C.png"), f4, width=8.5, height=8, dpi=300)

cat("wrote 4 result figures to", D, "\n")
