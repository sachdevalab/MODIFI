#!/usr/bin/env Rscript
# Fig 5b/5c regenerated on the revised strict-set high-conf linkages (rescued hosts).
# Reads the new mge_host_gc_cov_final.csv; same selection/test/style as count_linkages.R.
suppressMessages({library(readr); library(dplyr); library(ggplot2); library(tidyr)})
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr", repos="https://cloud.r-project.org")
suppressMessages(library(ggpubr))

O <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
gc_df <- read_csv(file.path(O, "mge_host_gc_cov_final.csv"), show_col_types = FALSE) %>%
  mutate(MGE_type = factor(MGE_type, levels = c("plasmid", "virus")))

cat("=== Cosine similarity ===\n")
print(gc_df %>% group_by(MGE_type) %>% summarise(mean=mean(cos_sim,na.rm=T), sd=sd(cos_sim,na.rm=T), n=n(), .groups="drop"))
print(gc_df %>% summarise(mean=mean(cos_sim,na.rm=T), sd=sd(cos_sim,na.rm=T), n=n()))

star_fun <- function(p) symnum(p, corr=FALSE, na=FALSE, cutpoints=c(0,0.001,0.01,0.05,0.1,1), symbols=c("***","**","*",".",""))
mytheme <- theme_bw() + theme(plot.title=element_text(hjust=.5,size=14), axis.title=element_text(size=12),
  axis.text=element_text(size=10), axis.text.x=element_text(angle=45,hjust=1), legend.position="top")

# ---- Panel b: linkages per sample by environment (>=4 samples) ----
env_counts <- gc_df %>% group_by(environment, sample, MGE_type, .drop=FALSE) %>% summarise(count=n(), .groups="drop") %>%
  group_by(environment) %>% filter(n_distinct(sample) >= 4) %>% ungroup()
env_order <- gc_df %>% count(environment, sort=TRUE) %>% pull(environment)
env_counts$environment <- factor(env_counts$environment, levels=env_order)
env_p <- env_counts %>% group_by(environment) %>% filter(length(unique(MGE_type))==2) %>%
  summarise(p=tryCatch(t.test(count~MGE_type)$p.value, error=function(e) NA_real_)) %>% mutate(star=star_fun(p))
env_p <- left_join(env_p, env_counts %>% group_by(environment) %>% summarise(y=max(count,na.rm=T)*1.1), by="environment")
cat("\n=== env p-values ===\n"); print(env_p)
p1 <- ggplot(env_counts, aes(environment, count, fill=MGE_type)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=.8), color="black", size=.5) +
  geom_jitter(position=position_jitterdodge(jitter.width=.2, dodge.width=.8), color="black", size=1.2, show.legend=FALSE) +
  scale_fill_discrete(name="MGE Type") + labs(x="Habitat", y="Linkages per Sample") + mytheme +
  geom_text(data=env_p, aes(environment, y, label=star), inherit.aes=FALSE, vjust=0, size=6)
ggsave(file.path(O,"fig5b_linkage_counts_env.pdf"), p1, width=3.5, height=6)
ggsave(file.path(O,"fig5b_linkage_counts_env.png"), p1, width=3.5, height=6, dpi=200)

# ---- Panel c: linkages per sample by host phylum (>=5 samples, top 7) ----
enough <- gc_df %>% group_by(host_phylum) %>% filter(n_distinct(sample) >= 5) %>% ungroup()
top7 <- enough %>% count(host_phylum, sort=TRUE) %>% head(7) %>% pull(host_phylum)
ph_counts <- enough %>% filter(host_phylum %in% top7) %>% mutate(host_phylum=gsub("^p__","",host_phylum)) %>%
  group_by(host_phylum, sample, MGE_type, .drop=FALSE) %>% summarise(count=n(), .groups="drop")
ph_counts$host_phylum <- factor(ph_counts$host_phylum, levels=gsub("^p__","",top7))
ph_p <- ph_counts %>% group_by(host_phylum) %>% filter(length(unique(MGE_type))==2) %>%
  summarise(p=tryCatch(t.test(count~MGE_type)$p.value, error=function(e) NA_real_)) %>% mutate(star=star_fun(p))
ph_p <- left_join(ph_p, ph_counts %>% group_by(host_phylum) %>% summarise(y=max(count,na.rm=T)*1.1), by="host_phylum")
cat("\n=== phylum p-values ===\n"); print(ph_p)
p2 <- ggplot(ph_counts, aes(host_phylum, count, fill=MGE_type)) +
  geom_boxplot(outlier.shape=NA, position=position_dodge(width=.8), color="black", size=.5) +
  geom_jitter(position=position_jitterdodge(jitter.width=.2, dodge.width=.8), color="black", size=1.2, show.legend=FALSE) +
  scale_fill_discrete(name="MGE Type") + labs(x="Host Phylum", y="Linkages per Sample") +
  mytheme + theme(legend.position="none") +
  geom_text(data=ph_p, aes(host_phylum, y, label=star), inherit.aes=FALSE, vjust=0, size=6)
ggsave(file.path(O,"fig5c_linkage_counts_phylum.pdf"), p2, width=7, height=5)
ggsave(file.path(O,"fig5c_linkage_counts_phylum.png"), p2, width=7, height=5, dpi=200)
cat("\nwrote fig5b_linkage_counts_env + fig5c_linkage_counts_phylum (pdf/png)\n")
