#!/usr/bin/env Rscript
# Figures for the all-linked-ECE analysis. Run with r_env Rscript.
suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(scales); library(patchwork)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
m <- read.csv(file.path(D,"all_linked_eces_annotated.tsv"), sep="\t", stringsAsFactors=FALSE)

env_lab <- c(infant_gut="infant gut", adult_gut="adult gut", mice_gut="mouse gut",
             cow_bioreactor="cow bioreactor", cow_rumen="cow rumen", soil="soil",
             sugarcane="sugarcane", mock="mock")
m$env <- factor(ifelse(m$environment %in% names(env_lab), env_lab[m$environment], m$environment),
                levels=env_lab)

## Fig A: known-vs-inferred host concordance
kn <- m %>% filter(known_plasmid=="True" | known_plasmid==TRUE)
top <- kn %>% count(inferred_host_genus, sort=TRUE) %>% head(8) %>% pull(inferred_host_genus)
coll <- function(x) ifelse(x %in% top & x!="", x, ifelse(x=="","unclassified","other"))
kn$ig <- coll(kn$inferred_host_genus); kn$rg <- coll(kn$ref_host_genus)
conc <- kn %>% count(ig, rg, name="n")
lv <- c(sort(unique(top)), "other","unclassified")
conc$ig <- factor(conc$ig, levels=rev(lv)); conc$rg <- factor(conc$rg, levels=lv)
conc$diag <- as.character(conc$ig)==as.character(conc$rg)
pA1 <- ggplot(conc, aes(rg, ig, fill=n)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_tile(data=subset(conc,diag), color="#d7301f", fill=NA, linewidth=1.1) +
  geom_text(aes(label=n), size=3) +
  scale_fill_gradient(low="#eef3f8", high="#08519c", name="plasmids") +
  labs(x="Documented host genus of closest known plasmid",
       y="Modification-inferred host genus",
       title="a  Inferred host vs known-plasmid host (n=84 known plasmids)",
       subtitle="red outline = agreement (on-diagonal); 74% agree at genus level") +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1,face="italic"),
        axis.text.y=element_text(face="italic"), panel.grid=element_blank(),
        plot.title=element_text(size=11), plot.subtitle=element_text(size=9))
agr <- kn %>% mutate(ok=ref_matches_inferred_genus=="yes") %>%
  group_by(env) %>% summarise(y=sum(ok), n=n(), .groups="drop") %>%
  filter(n>0) %>% mutate(rate=y/n, lab=paste0(y,"/",n))
pA2 <- ggplot(agr, aes(rate, env)) +
  geom_col(fill="#4292c6", width=0.7) +
  geom_text(aes(label=lab), hjust=-0.15, size=3) +
  scale_x_continuous(labels=percent, limits=c(0,1.12), breaks=c(0,.5,1)) +
  labs(x="genus-level agreement", y=NULL, title="b  Agreement by environment") +
  theme_classic(base_size=10) + theme(plot.title=element_text(size=11))
fA <- pA1 + pA2 + plot_layout(widths=c(1.6,1))
ggsave(file.path(D,"fig_linked_validation.pdf"), fA, width=11, height=5.2)
ggsave(file.path(D,"fig_linked_validation.png"), fA, width=11, height=5.2, dpi=300)

## Fig B: linked-ECE landscape across environments
pB1 <- ggplot(m, aes(env, fill=MGE_type)) +
  geom_bar(width=0.7) +
  scale_fill_manual(values=c(plasmid="#2171b5", virus="#238b45", novel="#bdbdbd"), name="ECE type") +
  labs(x=NULL, y="linked ECEs", title="a  Linked ECEs by environment") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1), plot.title=element_text(size=11))
pl <- m %>% filter(MGE_type=="plasmid")
pl$mob <- ifelse(grepl("conjugative", pl$predicted_mobility), "conjugative",
          ifelse(grepl("(^|;)mobilizable", pl$predicted_mobility), "mobilizable", "non-mobilizable"))
pl$mob <- factor(pl$mob, levels=c("conjugative","mobilizable","non-mobilizable"))
pB2 <- ggplot(pl, aes(env, fill=mob)) +
  geom_bar(position="fill", width=0.7) +
  scale_fill_manual(values=c(conjugative="#cb181d", mobilizable="#fd8d3c", `non-mobilizable`="#bdbdbd"),
                    name="Mobility") +
  scale_y_continuous(labels=percent) +
  labs(x=NULL, y="fraction of plasmids", title="b  Plasmid mobility by environment") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1), plot.title=element_text(size=11))
kf <- pl %>% mutate(known=known_plasmid=="True"|known_plasmid==TRUE) %>%
  group_by(env) %>% summarise(rate=mean(known), n=n(), .groups="drop") %>% filter(n>0)
pB3 <- ggplot(kf, aes(env, rate)) +
  geom_col(fill="#6a51a3", width=0.7) +
  geom_text(aes(label=paste0(round(100*rate),"%")), vjust=-0.4, size=3) +
  scale_y_continuous(labels=percent, expand=expansion(mult=c(0,0.12))) +
  labs(x=NULL, y="plasmids matching a known plasmid",
       title="c  Known-plasmid fraction by environment") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1), plot.title=element_text(size=11))
fB <- pB1 / pB2 / pB3
ggsave(file.path(D,"fig_linked_overview.pdf"), fB, width=7.5, height=10)
ggsave(file.path(D,"fig_linked_overview.png"), fB, width=7.5, height=10, dpi=300)
cat("wrote fig_linked_validation and fig_linked_overview to", D, "\n")
