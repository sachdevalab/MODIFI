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
# One label per column: REBASE homolog with its recognition motif appended in
# parentheses when there is one (no "n/a" for motif-less genes). HMM family, mod
# type and RM system are shown as annotation bars above.
homolog    <- sub("^\\[HMM\\] ", "", rownames(ann))          # REBASE homolog
has_motif  <- !(ann$motif == "" | is.na(ann$motif))
col_labels <- ifelse(has_motif, paste0(homolog, " (", ann$motif, ")"), homolog)
modtype    <- ifelse(ann$mod_type == "" | is.na(ann$mod_type), "unknown", ann$mod_type)
systype    <- ifelse(ann$system_type == "" | is.na(ann$system_type), "unknown", ann$system_type)
hmmfam     <- ann$hmm

mod_pal <- c(m6A = "#1b7837", m5C = "#762a83", m4C = "#e08214", unknown = "#dddddd")
sys_pal <- c("Type I" = "#4575b4", "Type II" = "#91bfdb", "Type IIG" = "#fee090",
             "Type III" = "#fc8d59", "Type IV" = "#d73027", unknown = "#dddddd")
hmm_levels <- unique(hmmfam)
hmm_pal <- structure(grDevices::hcl.colors(length(hmm_levels), "Dark 3"),
                     names = hmm_levels)

top_anno <- HeatmapAnnotation(
  `HMM family`  = hmmfam,
  `Mod type`    = modtype,
  `RM system`   = systype,
  col = list(`HMM family` = hmm_pal, `Mod type` = mod_pal, `RM system` = sys_pal),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "left",
  simple_anno_size = unit(3.5, "mm"),
  gap = unit(1, "mm"),
  show_legend = FALSE            # legends built manually and packed into one column
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
  cluster_columns = FALSE,                         # keep HMM-family column order
  column_order = colnames(pres),

  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 9),
  column_labels = col_labels,      # REBASE homolog (motif) on one line
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 90,

  top_annotation = top_anno,
  width  = unit(0.55 * ncol(pres), "cm"),
  height = unit(0.55 * nrow(pres), "cm"),

  show_heatmap_legend = FALSE     # legends built manually and packed into one column
)

# ---- all legends, stacked in a single column ----
lt <- gpar(fontsize = 9, fontface = "bold")
ll <- gpar(fontsize = 8)

gene_lgd <- Legend(title = "MTase gene", title_gp = lt, labels_gp = ll,
                   labels = c("absent", "present"),
                   legend_gp = gpar(fill = c("#f0f0f0", "#3690c0")))

mod_present <- intersect(names(mod_pal), unique(modtype))
mod_lgd <- Legend(title = "Mod type", title_gp = lt, labels_gp = ll,
                  labels = mod_present, legend_gp = gpar(fill = mod_pal[mod_present]))

sys_present <- intersect(names(sys_pal), unique(systype))
sys_lgd <- Legend(title = "RM system", title_gp = lt, labels_gp = ll,
                  labels = sys_present, legend_gp = gpar(fill = sys_pal[sys_present]))

hmm_lgd <- Legend(title = "HMM family", title_gp = lt, labels_gp = ll,
                  labels = hmm_levels, legend_gp = gpar(fill = hmm_pal[hmm_levels]))

act_lgd <- Legend(
  title = "Activity (motif)", title_gp = lt, labels_gp = ll,
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

packed <- packLegend(gene_lgd, act_lgd, mod_lgd, sys_lgd, hmm_lgd,
                     direction = "vertical", gap = unit(5, "mm"),
                     max_height = unit(30, "cm"))

w_in <- 0.55 * ncol(pres) / 2.54 + 5
h_in <- 0.55 * nrow(pres) / 2.54 + 5.5

for (dev in c("pdf", "png")) {
  f <- file.path(io, paste0("Kp248_1_mtase_presence_absence.", dev))
  if (dev == "pdf") pdf(f, width = w_in, height = h_in)
  else png(f, width = w_in, height = h_in, units = "in", res = 300)
  draw(ht, annotation_legend_side = "right", annotation_legend_list = list(packed))
  dev.off()
  cat(sprintf("wrote %s\n", f))
}
