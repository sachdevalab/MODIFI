#!/usr/bin/env Rscript
suppressMessages({library(ggplot2); library(dplyr); library(scales)})
D <- "/home/shuaiw/MODIFI/tmp/rev_figs/linked_eces"
s <- read.csv(file.path(D,"ece_validation_summary_strict.csv"), stringsAsFactors=FALSE)
catlev <- c("host-supported (agrees)","host mismatch",
            "known (strict), host not genus-resolved","no strict match")
pal <- c("host-supported (agrees)"="#238b45","host mismatch"="#d7301f",
         "known (strict), host not genus-resolved"="#9ecae1","no strict match"="#e0e0e0")
rate <- s %>% group_by(type) %>% summarise(
  total=sum(n),
  known=sum(n[category!="no strict match"]),
  resolved=sum(n[category %in% c("host-supported (agrees)","host mismatch")]),
  supp=sum(n[category=="host-supported (agrees)"]), .groups="drop")
rate$lab <- sprintf("%s\nknown %d/%d (%.0f%%);  host-supported %d/%d (%.0f%%)",
                    rate$type, rate$known, rate$total, 100*rate$known/rate$total,
                    rate$supp, rate$resolved, 100*rate$supp/rate$resolved)
lab_map <- setNames(rate$lab, rate$type)
s$category <- factor(s$category, levels=rev(catlev))
s$type <- factor(lab_map[s$type], levels=lab_map[c("Viruses (IMG/VR)","Plasmids (IMG/PR)")])

p <- ggplot(s, aes(x=n, y=type, fill=category)) +
  geom_col(width=0.62, color="white", linewidth=0.3) +
  geom_text(data=subset(s, n>=3), aes(label=n),
            position=position_stack(vjust=0.5), size=3.6, color="grey15") +
  scale_fill_manual(values=pal, breaks=catlev, name=NULL) +
  scale_x_continuous(expand=expansion(mult=c(0,0.03))) +
  labs(x="linked ECEs", y=NULL,
       title="Independent host validation of modification-linked ECEs",
       subtitle="Known = best hit to a metagenomic reference catalogue at >=90% identity & >=50% query coverage\n(IMG/VR for viruses, IMG/PR for plasmids); host-supported if the inferred genus matches ANY strict hit reference's genus-level host") +
  theme_classic(base_size=12) +
  theme(plot.title=element_text(size=13,face="bold"),
        plot.subtitle=element_text(size=9, color="grey30"),
        axis.text.y=element_text(size=10, lineheight=0.95),
        legend.position="bottom", legend.text=element_text(size=8.5)) +
  guides(fill=guide_legend(nrow=1, byrow=TRUE))
ggsave(file.path(D,"fig_ece_validation_summary.pdf"), p, width=12.5, height=4.8)
ggsave(file.path(D,"fig_ece_validation_summary.png"), p, width=12.5, height=4.8, dpi=300)
cat("wrote fig_ece_validation_summary (strict)\n")
