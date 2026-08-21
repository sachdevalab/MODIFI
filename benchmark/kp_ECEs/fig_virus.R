#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr); library(scales); library(patchwork)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
v <- read.csv(file.path(D,"virus_validation_categories.csv"), stringsAsFactors=FALSE)
v$category <- gsub("\\n"," ", v$category)

lev <- c("known virus, host confirmed","CRISPR-spacer host confirmed",
         "weak/distant hit (likely novel)","no reference hit (novel)")
cnt <- v %>% count(category) %>% mutate(category=factor(category, levels=rev(lev))) %>% arrange(category)
pal <- c("known virus, host confirmed"="#238b45","CRISPR-spacer host confirmed"="#66c2a4",
         "weak/distant hit (likely novel)"="#bdbdbd","no reference hit (novel)"="#e0e0e0")
pA <- ggplot(cnt, aes(n, category, fill=category)) +
  geom_col(width=0.72) + geom_text(aes(label=n), hjust=-0.3, size=3.4) +
  scale_fill_manual(values=pal, guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,0.12))) +
  labs(x="linked viruses (n=116)", y=NULL,
       title="a  Virus-linkage validation",
       subtitle="18/116 links independently checkable; 18/18 agree with the modification-inferred host") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(size=11,face="bold"), plot.subtitle=element_text(size=8.8))

env_lab <- c(cow_bioreactor="cow bioreactor", cow_rumen="cow rumen", infant_gut="infant gut",
             mice_gut="mouse gut", soil="soil")
cr <- read.csv("/home/shuaiw/MODIFI/benchmark/spacer/crispr_validation_breakdown.tsv", sep="\t")
cr <- cr %>% filter(type=="virus") %>%
  group_by(environment) %>% summarise(confirmed=sum(validated), total=n(), .groups="drop")
cr$env <- ifelse(cr$environment %in% names(env_lab), env_lab[cr$environment], cr$environment)
cr$env <- factor(cr$env, levels=cr$env[order(cr$confirmed)])
pB <- ggplot(cr, aes(confirmed, env)) +
  geom_col(fill="#66c2a4", width=0.7) +
  geom_text(aes(label=paste0(confirmed,"/",total)), hjust=-0.2, size=3.2) +
  scale_x_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x="CRISPR-spacer confirmed viruses", y=NULL, title="b  CRISPR confirmation by environment") +
  theme_classic(base_size=11) + theme(plot.title=element_text(size=11,face="bold"))

f <- pA / pB + plot_layout(heights=c(1,0.85))
ggsave(file.path(D,"fig_virus_validation.pdf"), f, width=8.5, height=6)
ggsave(file.path(D,"fig_virus_validation.png"), f, width=8.5, height=6, dpi=300)
cat("wrote fig_virus_validation\n")
