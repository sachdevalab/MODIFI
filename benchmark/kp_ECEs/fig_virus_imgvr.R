#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr); library(scales); library(patchwork)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
v <- read.csv(file.path(D,"virus_imgvr_validation.tsv"), sep="\t", stringsAsFactors=FALSE)
v <- v %>% filter(tier %in% c("strict","related"), ref_matches_inferred %in% c("yes","no"))

# genus concordance (collapse rare to "other")
top <- v %>% count(inferred_host_genus,sort=TRUE) %>% head(9) %>% pull(inferred_host_genus)
coll <- function(x) ifelse(x %in% top & x!="", x, "other")
v$ig <- coll(v$inferred_host_genus); v$rg <- coll(v$ref_host_genus)
conc <- v %>% count(ig,rg,name="n")
lv <- c(sort(setdiff(unique(c(v$ig,v$rg)),"other")),"other")
conc$ig <- factor(conc$ig, levels=rev(lv)); conc$rg <- factor(conc$rg, levels=lv)
conc$diag <- as.character(conc$ig)==as.character(conc$rg)
n_tot <- nrow(v); n_ok <- sum(v$ref_matches_inferred=="yes")
pA <- ggplot(conc, aes(rg, ig, fill=n)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_tile(data=subset(conc,diag), color="#d7301f", fill=NA, linewidth=1.1) +
  geom_text(aes(label=n), size=3) +
  scale_fill_gradient(low="#eef3f8", high="#08519c", name="viruses") +
  labs(x="IMG/VR predicted host genus", y="Modification-inferred host genus",
       title=sprintf("a  Virus host: inferred vs IMG/VR (n=%d genus-resolved)", n_tot),
       subtitle=sprintf("red outline = agreement; %d/%d (%.0f%%) agree at genus level",
                        n_ok, n_tot, 100*n_ok/n_tot)) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1,face="italic"),
        axis.text.y=element_text(face="italic"), panel.grid=element_blank(),
        plot.title=element_text(size=11), plot.subtitle=element_text(size=9))

# coverage leap: RefSeq vs IMG/VR
cov <- data.frame(
  db=factor(c("RefSeq viral","IMG/VR","RefSeq viral","IMG/VR"),levels=c("RefSeq viral","IMG/VR")),
  metric=c("viruses with hit","viruses with hit","genus-resolved host check","genus-resolved host check"),
  n=c(15,110,4,n_tot))
pB <- ggplot(cov, aes(n, db, fill=db)) +
  geom_col(width=0.7) + geom_text(aes(label=n), hjust=-0.3, size=3.2) +
  facet_wrap(~metric, ncol=1, scales="free_x") +
  scale_fill_manual(values=c("RefSeq viral"="#bdbdbd","IMG/VR"="#238b45"), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(x=NULL, y=NULL, title="b  IMG/VR vs RefSeq coverage (n=116 viruses)") +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(size=11), strip.text=element_text(size=9), strip.background=element_blank())

f <- pA + pB + plot_layout(widths=c(1.5,1))
ggsave(file.path(D,"fig_virus_imgvr_validation.pdf"), f, width=11, height=5.4)
ggsave(file.path(D,"fig_virus_imgvr_validation.png"), f, width=11, height=5.4, dpi=300)
cat("wrote fig_virus_imgvr_validation\n")
