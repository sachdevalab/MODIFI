#!/usr/bin/env Rscript
# Presence/absence heatmap of predicted methyltransferases (MTases) across the
# 15 Klebsiella pneumoniae 248_1 MAGs, with an activity overlay.
#
# Base cell:   gene present (dark) vs absent (light grey).
# Overlay dot: filled  = MTase recognition motif detected as modified in the MAG
#              open     = present but motif NOT detected (on/off - silent MTase)
#              (none)   = MTase has no REBASE recognition motif to assess
# Column labels = REBASE homolog (or HMM family); a text track above each column
# shows the recognition motif; two bars show modification type and RM system type.

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

io <- "/home/shuaiw/MODIFI/tmp/rev_figs/Klevsuekka_mtase"

# MAG (row) order transcribed from the published Fig 3e (top -> bottom), so this
# supplementary heatmap lines up row-for-row with the main-text figure.
FIG3E_ORDER <- c(
  "infant_27_13_C", "infant_1_1_C", "infant_25_32_C", "infant_8_9_C",
  "infant_19_5_L", "infant_20_6_C", "infant_3_4_C", "infant_2_6_C",
  "infant_14_14_C", "infant_3_25_C", "infant_2_11_C", "infant_4_8_C",
  "infant_15_6_C", "96plex_11_C", "infant_16_9_L"
)

pres <- as.matrix(read.csv(file.path(io, "mtase_presence_absence.csv"),
                           row.names = 1, check.names = FALSE))
act  <- as.matrix(read.csv(file.path(io, "mtase_activity.csv"),
                           row.names = 1, check.names = FALSE))
ann  <- read.csv(file.path(io, "mtase_column_annotation.csv"),
                 row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# align annotation to matrix columns
ann <- ann[colnames(pres), , drop = FALSE]

# ---- column labels ----
# Column keys are ortholog clusters; the cluster label already carries the
# representative REBASE homolog and its recognition motif, e.g. "M.Ksp2N3DamP (GATC)".
col_labels <- colnames(pres)
modtype    <- ifelse(ann$mod_type == "" | is.na(ann$mod_type), "unknown", ann$mod_type)
systype    <- ifelse(ann$system_type == "" | is.na(ann$system_type), "unknown", ann$system_type)

mod_pal <- c(m6A = "#1b7837", m5C = "#762a83", m4C = "#e08214", unknown = "#dddddd")
sys_pal <- c("Type I" = "#4575b4", "Type II" = "#91bfdb", "Type IIG" = "#fee090",
             "Type III" = "#fc8d59", "Type IV" = "#d73027", unknown = "#dddddd")

top_anno <- HeatmapAnnotation(
  `Mod type`    = modtype,
  `RM system`   = systype,
  col = list(`Mod type` = mod_pal, `RM system` = sys_pal),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "left",
  simple_anno_size = unit(3.5, "mm"),
  gap = unit(1, "mm"),
  show_legend = TRUE
)

# ---- activity overlay ----
overlay <- function(j, i, x, y, w, h, fill) {
  if (pres[i, j] == 1) {
    a <- act[i, j]
    if (!is.na(a) && a == 1) {                    # present & active -> filled dot
      grid.circle(x, y, r = unit(1.4, "mm"),
                  gp = gpar(fill = "black", col = NA))
    } else if (!is.na(a) && a == 0) {             # present but silent -> open dot
      grid.circle(x, y, r = unit(1.4, "mm"),
                  gp = gpar(fill = NA, col = "black", lwd = 1))
    }                                             # NA (no motif) -> no dot
  }
}

col_fun <- c("0" = "#f0f0f0", "1" = "#3690c0")

ht <- Heatmap(
  pres,
  name = "MTase gene",
  col = col_fun,
  rect_gp = gpar(col = "white", lwd = 1),
  cell_fun = overlay,

  cluster_rows = FALSE,                            # fixed MAG order (Fig 3e)
  row_order = FIG3E_ORDER,
  cluster_columns = TRUE,
  clustering_distance_columns = "binary",
  clustering_method_columns = "average",

  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_labels = col_labels,
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 90,

  top_annotation = top_anno,
  width  = unit(0.55 * ncol(pres), "cm"),
  height = unit(0.55 * nrow(pres), "cm"),

  heatmap_legend_param = list(
    title = "MTase gene",
    at = c(0, 1), labels = c("absent", "present"),
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

# manual legend for the activity overlay
act_lgd <- Legend(
  title = "Activity (motif)",
  title_gp = gpar(fontsize = 9, fontface = "bold"),
  labels_gp = gpar(fontsize = 8),
  labels = c("detected", "not detected", "no REBASE motif"),
  graphics = list(
    function(x, y, w, h) grid.circle(x, y, r = unit(1.4, "mm"),
                                     gp = gpar(fill = "black", col = NA)),
    function(x, y, w, h) grid.circle(x, y, r = unit(1.4, "mm"),
                                     gp = gpar(fill = NA, col = "black", lwd = 1)),
    function(x, y, w, h) grid.rect(x, y, unit(2.8, "mm"), unit(2.8, "mm"),
                                   gp = gpar(fill = "#3690c0", col = "white"))
  )
)

w_in <- 0.55 * ncol(pres) / 2.54 + 6
h_in <- 0.55 * nrow(pres) / 2.54 + 4

for (dev in c("pdf", "png")) {
  f <- file.path(io, paste0("Kp248_1_mtase_presence_absence.", dev))
  if (dev == "pdf") pdf(f, width = w_in, height = h_in)
  else png(f, width = w_in, height = h_in, units = "in", res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
       annotation_legend_list = list(act_lgd))
  dev.off()
  cat(sprintf("wrote %s\n", f))
}
