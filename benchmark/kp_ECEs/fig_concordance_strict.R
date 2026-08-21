#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr); library(scales); library(patchwork)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"

concordance <- function(tsv, refname, hi, hicol, fillhi, covrows, outbase){
  p <- read.csv(file.path(D,tsv), sep="\t", stringsAsFactors=FALSE)
  p <- p %>% filter(category %in% c("host-supported (agrees)","host mismatch"),
                    ref_host_genus!="", inferred_host_genus!="")
  top <- p %>% count(inferred_host_genus,sort=TRUE) %>% head(9) %>% pull(inferred_host_genus)
  coll <- function(x) ifelse(x %in% top & x!="", x, "other")
  p$ig <- coll(p$inferred_host_genus); p$rg <- coll(p$ref_host_genus)
  conc <- p %>% count(ig,rg,name="n")
  lv <- c(sort(setdiff(unique(c(p$ig,p$rg)),"other")),"other")
  conc$ig <- factor(conc$ig, levels=rev(lv)); conc$rg <- factor(conc$rg, levels=lv)
  conc$diag <- as.character(conc$ig)==as.character(conc$rg)
  n_tot <- nrow(p); n_ok <- sum(p$category=="host-supported (agrees)")
  pA <- ggplot(conc, aes(rg, ig, fill=n)) +
    geom_tile(color="white", linewidth=0.6) +
    geom_tile(data=subset(conc,diag), color="#d7301f", fill=NA, linewidth=1.1) +
    geom_text(aes(label=n), size=2.9) +
    scale_fill_gradient(low="#eef3f8", high=fillhi, name=hi) +
    labs(x=sprintf("%s predicted host genus",refname), y="Modification-inferred host genus",
         title=sprintf("a  %s host: inferred vs %s (strict, n=%d)", hi, refname, n_tot),
         subtitle=sprintf("known = >=90%% id & >=50%% cov; %d/%d (%.0f%%) agree at genus level",
                          n_ok,n_tot,100*n_ok/n_tot)) +
    theme_minimal(base_size=10) +
    theme(axis.text.x=element_text(angle=35,hjust=1,face="italic"),
          axis.text.y=element_text(face="italic"), panel.grid=element_blank(),
          plot.title=element_text(size=11), plot.subtitle=element_text(size=8.8))
  cov <- covrows
  pB <- ggplot(cov, aes(n, db, fill=db)) +
    geom_col(width=0.7) + geom_text(aes(label=n), hjust=-0.3, size=3.2) +
    facet_wrap(~metric, ncol=1, scales="free_x") +
    scale_fill_manual(values=setNames(c("#bdbdbd",fillhi),levels(cov$db)), guide="none") +
    scale_x_continuous(expand=expansion(mult=c(0,0.2))) +
    labs(x=NULL,y=NULL,title="b  Coverage vs clinical/RefSeq DB") +
    theme_classic(base_size=10) +
    theme(plot.title=element_text(size=10.5), strip.text=element_text(size=9), strip.background=element_blank())
  f <- pA + pB + plot_layout(widths=c(1.5,1))
  ggsave(file.path(D,paste0(outbase,".pdf")), f, width=11.5, height=6)
  ggsave(file.path(D,paste0(outbase,".png")), f, width=11.5, height=6, dpi=300)
}

vcov <- data.frame(db=factor(c("RefSeq viral","IMG/VR","RefSeq viral","IMG/VR"),levels=c("RefSeq viral","IMG/VR")),
                   metric=c("known (strict)","known (strict)","genus-resolved host","genus-resolved host"),
                   n=c(4,42,4,30))
concordance("virus_imgvr_strict_validation.tsv","IMG/VR","viruses","Host taxonomy","#238b45",vcov,"fig_virus_imgvr_validation")
pcov <- data.frame(db=factor(c("mob_suite/RefSeq","IMG/PR","mob_suite/RefSeq","IMG/PR"),levels=c("mob_suite/RefSeq","IMG/PR")),
                   metric=c("known (strict)","known (strict)","genus-resolved host","genus-resolved host"),
                   n=c(84,135,62,114))
concordance("plasmid_imgpr_strict_validation.tsv","IMG/PR","plasmids","host_taxonomy","#6a51a3",pcov,"fig_plasmid_imgpr_validation")
cat("re-rendered strict concordance figures\n")
