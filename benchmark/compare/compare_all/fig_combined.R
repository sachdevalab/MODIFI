#!/usr/bin/env Rscript
# Combined figure: (top) per-contig motif-count correlation for each tool pair;
# (bottom) non-redundant motif count per tool + pairwise containment Jaccard.
# No panel titles; panels tagged a-e by patchwork.
suppressMessages({library(ggplot2); library(patchwork)})
FIG <- "/home/shuaiw/MODIFI/tmp/rev_figs/compare_all_meta"
COL <- c("MODIFI"="#D55E00","ipdSummary"="#E69F00","fibertools"="#009E73")
lev <- c("MODIFI","ipdSummary","fibertools")

## ---- top: correlation scatters ----
d <- read.csv(file.path(FIG,"fig_percontig_counts.csv"), stringsAsFactors=FALSE)
lim <- c(0, max(d$MODIFI,d$ipdSummary,d$fibertools)+2)
scat <- function(xc,yc){
  r <- cor(d[[xc]], d[[yc]])
  ggplot(d, aes(.data[[xc]], .data[[yc]])) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
    geom_point(position=position_jitter(0.15,0.15), alpha=0.45, size=1.5, color="#2171b5") +
    annotate("text", x=lim[1]+0.02*lim[2], y=lim[2],
             label=sprintf("r = %.2f", r), hjust=0, vjust=1, size=3.8) +
    coord_equal(xlim=lim, ylim=lim) +
    labs(x=xc, y=yc) +
    theme_classic(base_size=12)
}

## ---- bottom-left: non-redundant motif count ----
cnt <- read.csv(file.path(FIG,"fig_bars.counts.csv"), stringsAsFactors=FALSE)
cnt$tool <- sub("MODIFI-sub","MODIFI",cnt$tool)
cnt$tool <- factor(cnt$tool, levels=lev)
pcount <- ggplot(cnt, aes(tool, n_nonredundant, fill=tool)) +
  geom_col(width=0.68, color="grey20", linewidth=0.3) +
  geom_text(aes(label=n_nonredundant), vjust=-0.4, size=4, fontface="bold") +
  scale_fill_manual(values=COL, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.13))) +
  labs(x=NULL, y="Non-redundant motifs") +
  theme_classic(base_size=12) +
  theme(axis.text.x=element_text(angle=20, hjust=1))

## ---- bottom-right: pairwise Jaccard on dereplicated motif sets ----
j <- read.csv(file.path(FIG,"fig_jaccard_drep.csv"), stringsAsFactors=FALSE)
j$i1 <- match(j$tool1,lev); j$i2 <- match(j$tool2,lev)
j <- j[j$i1 < j$i2, ]
j$pair <- paste(j$tool1, j$tool2, sep="  vs\n")
j$pair <- factor(j$pair, levels=j$pair[order(-j$J_dereplicated)])
pjac <- ggplot(j, aes(pair, J_dereplicated)) +
  geom_col(width=0.6, fill="#3182bd", color="grey20", linewidth=0.3) +
  geom_text(aes(label=sprintf("%.2f", J_dereplicated)), vjust=-0.35, size=4) +
  scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,0.25), expand=expansion(mult=c(0,0.02))) +
  labs(x=NULL, y="Jaccard") +
  theme_classic(base_size=12)

top <- scat("MODIFI","ipdSummary") | scat("MODIFI","fibertools") | scat("ipdSummary","fibertools")
bottom <- pcount | pjac
p <- (top / bottom) + plot_layout(heights=c(1,1)) + plot_annotation(tag_levels="a") &
  theme(plot.tag = element_text(face="bold", size=15))

ggsave(file.path(FIG,"fig_motif_combined.png"), p, width=11, height=8, dpi=200)
ggsave(file.path(FIG,"fig_motif_combined.pdf"), p, width=11, height=8)
cat("wrote fig_motif_combined.png/.pdf\n")
