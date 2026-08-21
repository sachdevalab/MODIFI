#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr); library(scales); library(patchwork)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
p <- read.csv(file.path(D,"plasmid_imgpr_validation.tsv"), sep="\t", stringsAsFactors=FALSE)
p <- p %>% filter(tier %in% c("strict","related"), ref_matches_inferred %in% c("yes","no"))

top <- p %>% count(inferred_host_genus,sort=TRUE) %>% head(9) %>% pull(inferred_host_genus)
coll <- function(x) ifelse(x %in% top & x!="", x, "other")
p$ig <- coll(p$inferred_host_genus); p$rg <- coll(p$ref_host_genus)
conc <- p %>% count(ig,rg,name="n")
lv <- c(sort(setdiff(unique(c(p$ig,p$rg)),"other")),"other")
conc$ig <- factor(conc$ig, levels=rev(lv)); conc$rg <- factor(conc$rg, levels=lv)
conc$diag <- as.character(conc$ig)==as.character(conc$rg)
n_tot <- nrow(p); n_ok <- sum(p$ref_matches_inferred=="yes")
pA <- ggplot(conc, aes(rg, ig, fill=n)) +
  geom_tile(color="white", linewidth=0.6) +
  geom_tile(data=subset(conc,diag), color="#d7301f", fill=NA, linewidth=1.1) +
  geom_text(aes(label=n), size=2.9) +
  scale_fill_gradient(low="#eef3f8", high="#6a51a3", name="plasmids") +
  labs(x="IMG/PR predicted host genus", y="Modification-inferred host genus",
       title=sprintf("a  Plasmid host: inferred vs IMG/PR (n=%d genus-resolved)", n_tot),
       subtitle=sprintf("%d/%d (%.0f%%) agree at genus level; off-diagonal reflects plasmid host promiscuity",
                        n_ok, n_tot, 100*n_ok/n_tot)) +
  theme_minimal(base_size=10) +
  theme(axis.text.x=element_text(angle=35,hjust=1,face="italic"),
        axis.text.y=element_text(face="italic"), panel.grid=element_blank(),
        plot.title=element_text(size=11), plot.subtitle=element_text(size=8.6))
cov <- data.frame(
  db=factor(rep(c("mob_suite / RefSeq","IMG/PR"),2),levels=c("mob_suite / RefSeq","IMG/PR")),
  metric=c("plasmids with known match","plasmids with known match","genus-resolved host check","genus-resolved host check"),
  n=c(84,254,84,n_tot))
pB <- ggplot(cov, aes(n, db, fill=db)) +
  geom_col(width=0.7) + geom_text(aes(label=n), hjust=-0.3, size=3.2) +
  facet_wrap(~metric, ncol=1, scales="free_x") +
  scale_fill_manual(values=c("mob_suite / RefSeq"="#bdbdbd","IMG/PR"="#6a51a3"), guide="none") +
  scale_x_continuous(expand=expansion(mult=c(0,0.2))) +
  labs(x=NULL,y=NULL,title="b  IMG/PR vs clinical-DB coverage (n=263 plasmids)") +
  theme_classic(base_size=10) +
  theme(plot.title=element_text(size=10.5), strip.text=element_text(size=9), strip.background=element_blank())
f <- pA + pB + plot_layout(widths=c(1.5,1))
ggsave(file.path(D,"fig_plasmid_imgpr_validation.pdf"), f, width=11.5, height=6)
ggsave(file.path(D,"fig_plasmid_imgpr_validation.png"), f, width=11.5, height=6, dpi=300)
cat("wrote fig_plasmid_imgpr_validation\n")
