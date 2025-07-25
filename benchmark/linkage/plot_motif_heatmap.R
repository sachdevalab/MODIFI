## load the csv file as a matrix
library(ggplot2)
library(aplot)


pdf(file="motif_heatmap.pdf", width=12, height=10, onefile=FALSE)


df<-read.table("heatmap.csv", sep=",", header=TRUE)
df$value <- as.numeric(as.character(df$value))


p1<-ggplot(data = df, aes(x=motif, y=seq, fill=value)) + 
  geom_tile()+
   theme(axis.text.x = element_text(size = 8))+
   scale_x_discrete(guide = guide_axis(angle = 90))+
 geom_tile(color = "white")+
 scale_fill_gradient2(
    space = "Lab",   #midpoint = 0.3, limit = c(0,0.5),
    name="Freq") +
  theme_classic()

    # theme(
    # axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    # panel.grid.major = element_blank(),
    # panel.border = element_blank(),
    # panel.background = element_blank(),
    # axis.ticks = element_blank())

p1

# df2<-read.csv("../hla/hla_color.csv",header=T, sep=",")
# df2$y<-factor(df2$y,levels = rev(df2$y))
# p2<-ggplot(df2,aes(x=x,y=y))+
#   geom_tile(aes(fill=group))+
#   scale_x_continuous(expand = c(0,0))+
#   theme(panel.background = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         legend.position = "left",
#         legend.title = element_blank())+
#   scale_fill_manual(values = c("#d9e6eb", "#9fc3d5", "#8f96bd", "#2a347a", "#d6d69b"))
# #         

# p <- p1%>%
#   insert_left(p2,width = 0.02) %>%
# p

dev.off()