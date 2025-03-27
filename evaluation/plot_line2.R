library(ggplot2)


# [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
# [8] "#666666"

pdf(file="borg/customized/line_plot3.pdf", width=16, height=12, onefile=FALSE)

df<-read.table("borg/customized/bin_df3.csv", sep=",", header=TRUE)

head (df)

p <- ggplot(df, aes(x=pos, y=IPD)) +
  geom_line(size=0.5) +
  geom_point(aes(size = score, color = feature_type)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


vysg<-p+facet_wrap(~strand, scales = "free", nrow = 8)+
xlab('')+
ylab('IPD')+
  theme(legend.position = "none")

vysg

dev.off()