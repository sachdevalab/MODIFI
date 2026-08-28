#!/usr/bin/env Rscript
# Independent IMG-catalogue host validation of the FINAL strict linkage set (317 linked ECEs).
# Copy of benchmark/kp_ECEs/fig_summary.R; input = new summary CSV; KNOWN = cluster_MGE.py ANI>=95 & qcov>=85.
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/ece_anno"
s <- read.csv(file.path(D,"ece_validation_summary_strict.csv"), stringsAsFactors=FALSE)
catlev <- c("host-supported (agrees)","host mismatch","no comparable reference")
pal <- c("host-supported (agrees)"="#238b45","host mismatch"="#d7301f","no comparable reference"="#e0e0e0")
rate <- s %>% group_by(type) %>% summarise(
  total=sum(n),
  known=sum(n[category!="no strict match"]),
  resolved=sum(n[category %in% c("host-supported (agrees)","host mismatch")]),
  supp=sum(n[category=="host-supported (agrees)"]), .groups="drop")
rate$lab <- rate$type
lab_map <- setNames(rate$lab, rate$type)
# grey = no strict IMG hit OR a strict hit whose reference lacks a genus-level host
s$category[s$category %in% c("known (strict), host not genus-resolved","no strict match")] <- "no comparable reference"
s <- s %>% group_by(type, category) %>% summarise(n=sum(n), .groups="drop")
s$category <- factor(s$category, levels=rev(catlev))
s$type <- factor(lab_map[s$type], levels=lab_map[c("Viruses (IMG/VR)","Plasmids (IMG/PR)")])

p <- ggplot(s, aes(x=n, y=type, fill=category)) +
  geom_col(width=0.62, color="white", linewidth=0.3) +
  geom_text(data=subset(s, n>=3), aes(label=n),
            position=position_stack(vjust=0.5), size=3.6, color="grey15") +
  scale_fill_manual(values=pal, breaks=catlev, name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,0.03))) +
  labs(x="linked ECEs", y=NULL,
       subtitle="Each linked ECE searched against IMG/VR (viruses) or IMG/PR (plasmids); a reference match requires ANI>=95% & >=85% query coverage (same stringency as ECE clustering).\nHost-supported = inferred host genus matches the genus-level host of any matched reference; host mismatch = a matched reference has a genus-level host but it disagrees.") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(size=13,face="bold"),
        plot.subtitle=element_text(size=9, color="grey30"),
        axis.text.y=element_text(size=10, lineheight=0.95),
        legend.position="bottom", legend.text=element_text(size=8.5)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(file.path(D,"img_host_validation_final.pdf"), p, width=12.5, height=4.8)
ggsave(file.path(D,"img_host_validation_final.png"), p, width=12.5, height=4.8, dpi=300)
cat("wrote img_host_validation_final (strict, ANI>=95 & qcov>=85)\n")
