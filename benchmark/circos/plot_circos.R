library(circlize)

# Simulate a circular bacterial genome of 5 Mbp
genome_length <- 94795
num_motifs <- 30
pdf("../../tmp/circos_plot.pdf")

# Generate random methylated motif positions and types
# set.seed(42)
# motif_positions <- sample(1:genome_length, num_motifs)
# motif_types <- sample(c("m6A", "m5C"), num_motifs, replace = TRUE)
# motif_colors <- ifelse(motif_types == "m6A", "red", "blue")

# # Create data frame for motifs
# motif_df <- data.frame(
#   chr = "chr1",
#   start = motif_positions - 5000,
#   end = motif_positions + 5000,
#   color = motif_colors,
#   label = motif_types
# )
motif_df <- read.csv("./motif_df.csv")

head(motif_df)

# Define genome
genome_df <- data.frame(chr = "E_coli_H10407_2", start = 0, end = genome_length)

# Initialize circos
circos.clear()
circos.par(start.degree = 0, gap.degree = 0, track.margin = c(0.01, 0.01))
circos.initialize(factors = genome_df$chr, xlim = genome_df[, c("start", "end")])

# Draw a visible black genome ring
circos.trackPlotRegion(
  ylim = c(0, 1),
  track.height = 0.05,
  panel.fun = function(x, y) {
    circos.rect(CELL_META$xlim[1], 0, CELL_META$xlim[2], 1, col = "grey90", border = "black")
  },
  bg.border = NA
)

# # Overlay methylated motifs as colored rectangles
# for (i in 1:nrow(motif_df)) {
#   circos.rect(
#     xleft = motif_df$start[i], xright = motif_df$end[i],
#     ybottom = 0, ytop = 1,
#     col = motif_df$color[i], border = NA,
#     sector.index = motif_df$chr[i]
#   )
#   # circos.text(
#   #   x = (motif_df$start[i] + motif_df$end[i]) / 2,
#   #   y = 1.2,
#   #   labels = motif_df$label[i],
#   #   col = motif_df$color[i],
#   #   cex = 0.6,
#   #   facing = "clockwise",
#   #   niceFacing = TRUE,
#   #   sector.index = motif_df$chr[i]
#   # )
# }



## extract a new df with the motif_df has motif GATC
motif1_df <- motif_df[motif_df$motif == "GATC", ] 
motif2_df <- motif_df[motif_df$motif == "CTTCAG", ] 
motif3_df <- motif_df[motif_df$motif == "AGCANNNNNNCCT", ] 
motif4_df <- motif_df[motif_df$motif == "CAAYNNNNNCTGC", ] 

# Add a separate track for motifs
circos.genomicTrack(
  motif1_df,
  ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    i <- getI(...)
    circos.rect(
      xleft = region$start,
      xright = region$end,
      ybottom = 0,
      ytop = 1,
      col = motif_df$color,
      border = NA
    )
  },
  track.height = 0.05,
  bg.border = NA
)


# Add a separate track for motifs
circos.genomicTrack(
  motif2_df,
  ylim = c(0, 1),
  panel.fun = function(region, value, ...) {
    i <- getI(...)
    circos.rect(
      xleft = region$start,
      xright = region$end,
      ybottom = 0,
      ytop = 1,
      col = motif_df$color,
      border = NA
    )
  },
  track.height = 0.05,
  bg.border = NA
)




title("Simulated Bacterial Genome with Methylated Motifs")
dev.off()


