#!/usr/bin/env Rscript
# Annotated version of the Methanoperedens / Borg modification heatmap (Figure S14).
#
# Adds a left-side annotation with three color-keyed bars, per reviewer comment #9:
#   1. Type              - element type (Mp host / Borg / mini_Borg / Mini_Chr / Mp_Virus)
#   2. Associated host   - host group derived from the motif-profile NETWORK (same edges and
#                          cutoff as profile_borg.py: Pearson corr > 0.7 AND >=1 shared motif),
#                          so a host and the elements linked to it in the network share a color.
#   3. Sample            - source sample (replaces the old per-sample coloring of label text)
#
# Row-label text is plain black. The heatmap value legend is shown ("Modification fraction").
#
# This is a NEW script and writes NEW outputs; it does not modify plot_borg_heatmap.R or
# overwrite borg_profile_heatmap.pdf.
#
# Run with the r_env conda R:
#   /home/shuaiw/miniconda3/envs/r_env/bin/Rscript plot_borg_heatmap_annotated.R

library(tidyr)
library(dplyr)
library(grid)
library(ComplexHeatmap)
suppressPackageStartupMessages(library(circlize))
library(RColorBrewer)
library(igraph)

# Build a plotmath expression vector for legend labels: italicize the "Methanoperedens[ sp...]"
# binomial while leaving tags like "<GTDB>" upright. Non-Methanoperedens labels are wrapped as
# plain string constants so underscores are shown literally (not turned into subscripts).
italicize_mp <- function(labels) {
  exprs <- vector("expression", length(labels))
  for (i in seq_along(labels)) {
    lab <- labels[i]
    if (grepl("^Methanoperedens", lab)) {
      m <- regmatches(lab, regexec("^(Methanoperedens[^<]*?)\\s*(<.*>)?$", lab))[[1]]
      genus_part <- trimws(m[2])
      tail_part  <- if (length(m) >= 3 && nzchar(m[3])) trimws(m[3]) else ""
      if (nzchar(tail_part)) {
        exprs[[i]] <- bquote(italic(.(genus_part)) * " " * .(tail_part))
      } else {
        exprs[[i]] <- bquote(italic(.(genus_part)))
      }
    } else {
      exprs[[i]] <- bquote(.(lab))  # plain string, rendered literally
    }
  }
  exprs
}

# ---- IUPAC-aware motif deduplication helpers ----
.IUPAC <- list(A="A", C="C", G="G", T="T", R="AG", Y="CT", S="GC", W="AT", K="GT",
               M="AC", B="CGT", D="AGT", H="ACT", V="ACG", N="ACGT")
.iset <- function(ch) strsplit(.IUPAC[[ch]], "")[[1]]
.rc_iupac <- function(s) {
  comp <- c(A="T", T="A", C="G", G="C", R="Y", Y="R", S="S", W="W", K="M", M="K",
            B="V", V="B", D="H", H="D", N="N")
  ch <- strsplit(s, "")[[1]]
  paste(rev(comp[ch]), collapse = "")
}
# offsets of each base relative to the modified position (0-based)
.motif_offsets <- function(bare, modpos) {
  ch <- strsplit(bare, "")[[1]]
  setNames(lapply(seq_along(ch), function(i) .iset(ch[i])),
           as.character(seq_along(ch) - 1 - modpos))
}
# TRUE if motif A generalizes motif B (every site matching B also matches A)
.covers <- function(A, B) {
  full <- c("A", "C", "G", "T")
  for (o in names(A)) {
    sb <- if (!is.null(B[[o]])) B[[o]] else full
    if (!all(sb %in% A[[o]])) return(FALSE)
  }
  TRUE
}
# Motif dedup by redundancy (containment / reverse-complement), always keeping the motif
# with MORE modified sites. For any pair where one motif is a subset of the other (in either
# orientation), drop the one with fewer modified sites (`weight`); ties break by specificity
# then name. This collapses redundant families to their best-supported representative and so
# retains high-signal degenerate motifs (e.g. YCTB, ACC) rather than dropping them.
dedup_motifs_by_support <- function(motif_ids, weight) {
  parse1 <- function(m) {
    p <- strsplit(m, "_")[[1]]
    bare <- p[1]; modpos <- as.integer(p[2])
    w <- as.numeric(weight[m]); if (is.na(w)) w <- 0
    list(off = .motif_offsets(bare, modpos),
         offrc = .motif_offsets(.rc_iupac(bare), nchar(bare) - 1 - modpos),
         spec = sum(strsplit(bare, "")[[1]] != "N"), w = w)
  }
  I <- lapply(motif_ids, parse1); names(I) <- motif_ids
  related <- function(a, b) .covers(I[[a]]$off, I[[b]]$off) ||
                            .covers(I[[b]]$off, I[[a]]$off) ||
                            .covers(I[[a]]$off, I[[b]]$offrc) ||
                            .covers(I[[b]]$off, I[[a]]$offrc)
  loses <- function(a, b) {  # a is dropped in favor of b
    if (I[[a]]$w   != I[[b]]$w)   return(I[[a]]$w   < I[[b]]$w)      # fewer modified sites
    if (I[[a]]$spec != I[[b]]$spec) return(I[[a]]$spec < I[[b]]$spec) # then less specific
    a > b
  }
  drop <- character(0)
  for (a in motif_ids) {
    for (b in motif_ids) {
      if (a == b) next
      if (related(a, b) && loses(a, b)) { drop <- c(drop, a); break }
    }
  }
  setdiff(motif_ids, unique(drop))
}

# Concise per-genome labels for the transposed (single-sample) panel, so the long contig
# names do not dominate the axis. Borg/mini_Borg keep their named element; Mini_Chr keeps
# its size; hosts get "Mp <id>".
short_genome_label <- function(contig, type, borg_ref) {
  if (type %in% c("Borg", "mini_Borg")) {
    lab <- borg_ref
    lab <- sub("[-_](SR-VP|SRVP)[0-9]*_.*$", "", lab)  # drop "_SR-VP.../-SR-VP..." sample suffix
    lab <- sub("_scaf(fold)?.*$", "", lab)              # drop "_scaffold.../_scaf..." suffix
    if (!nzchar(lab) || is.na(lab)) lab <- borg_ref
    return(lab)
  }
  if (type == "Mini_Chr") {
    m <- regmatches(contig, regexpr("[0-9]+(-[0-9]+)?(kb|MB|Mb)", contig))
    if (length(m) > 0 && nzchar(m)) return(paste0("Mini_Chr_", m))
    return("Mini_Chr")
  }
  # Mp / Mp_Virus: trailing numeric contig id
  m <- regmatches(contig, regexpr("[0-9]+_[LC]$", contig))
  id <- if (length(m) > 0 && nzchar(m)) m else contig
  paste0(type, "_", id)
}

personal_plot_annotated <- function(profile_df, plot_name,
                                    network_edge_file,
                                    network_colors_file,
                                    sample_name_file = NULL,
                                    plot_height = 11,
                                    motif_min_fraction = 0,
                                    motif_min_sites = 0,
                                    motif_min_sites_strong = Inf,
                                    site_count_file = NULL,
                                    show_motif_names = FALSE,
                                    dedup_motifs = FALSE,
                                    transpose = FALSE,
                                    plot_width = 17) {

  # ---- Build contig x motif matrix (same as plot_borg_heatmap.R) ----
  pivot_df <- profile_df %>%
    select(contig, motifString, fraction) %>%
    pivot_wider(names_from = motifString, values_from = fraction, values_fill = 0)

  contig_names <- pivot_df$contig

  # Per-contig type (Genome) and BORG_Ref, captured before clustering
  contig_meta <- profile_df %>%
    group_by(contig) %>%
    summarise(Genome = first(Genome), BORG_Ref = first(BORG_Ref), .groups = 'drop')
  contig_to_genome <- setNames(contig_meta$Genome, contig_meta$contig)
  contig_to_borg   <- setNames(contig_meta$BORG_Ref, contig_meta$contig)

  pivot_matrix <- as.matrix(pivot_df %>% select(-contig))
  rownames(pivot_matrix) <- contig_names

  # Row labels: plain black, concise. Borg / mini_Borg keep their named element (BORG_Ref);
  # hosts get "Mp_<id>", Mini_Chr keeps its size, via short_genome_label().
  y_labels_all <- mapply(short_genome_label, contig_names,
                         contig_to_genome[contig_names], contig_to_borg[contig_names])

  # ---- Sample of origin per contig ----
  extract_sample_name <- function(contig_name) {
    parts <- strsplit(contig_name, "_")[[1]]
    if (length(parts) >= 8) paste(parts[1:8], collapse = "_") else contig_name
  }
  sample_name_map <- NULL
  if (!is.null(sample_name_file) && file.exists(sample_name_file)) {
    sndf <- read.csv(sample_name_file, stringsAsFactors = FALSE)
    if (all(c("contig", "sample") %in% colnames(sndf))) {
      sample_name_map <- sndf %>% distinct(contig, sample) %>% tibble::deframe()
    }
  }
  if (!is.null(sample_name_map)) {
    sample_names <- sapply(contig_names, function(c) {
      v <- sample_name_map[c]
      if (is.na(v) || length(v) == 0) extract_sample_name(c) else v
    })
  } else {
    sample_names <- sapply(contig_names, extract_sample_name)
  }
  names(sample_names) <- contig_names

  # ---- Type (Genome) per contig ----
  type_names <- contig_to_genome[contig_names]
  names(type_names) <- contig_names

  # ---- Linkage groups from the NETWORK (same cutoff as profile_borg.py) ----
  # Edges already encode Pearson corr > 0.7 AND >=1 shared motif. Connected components
  # (restricted to the heatmap contigs) are the sets of genomes thought to be linked
  # together. We keep the component id per contig (used to split the heatmap into blocks)
  # and a coarser label for the color bar:
  #   - a component that contains an Mp host  -> that Methanoperedens lineage ("host")
  #   - a component of only Borg/ECE elements -> "Borg or ECE-linked"        ("borg_ece")
  #   - a singleton                           -> "Unlinked"                  ("unlinked")
  edge_df <- read.csv(network_edge_file, stringsAsFactors = FALSE)
  edge_df <- edge_df %>%
    filter(Contig1 %in% contig_names, Contig2 %in% contig_names)

  g <- graph_from_data_frame(
    d = edge_df[, c("Contig1", "Contig2")],
    vertices = data.frame(name = contig_names, stringsAsFactors = FALSE),
    directed = FALSE
  )
  comp <- components(g)
  comp_id <- comp$membership[contig_names]
  names(comp_id) <- contig_names

  linkage_label <- rep("Unlinked", length(contig_names))
  linkage_class <- rep("unlinked", length(contig_names))
  names(linkage_label) <- contig_names
  names(linkage_class) <- contig_names
  for (cid in unique(comp_id)) {
    members <- contig_names[comp_id == cid]
    if (length(members) < 2) next  # singleton -> Unlinked
    mp_members <- members[type_names[members] == "Mp"]
    if (length(mp_members) > 0) {
      lineage <- names(sort(table(contig_to_borg[mp_members]), decreasing = TRUE))[1]
      linkage_label[members] <- lineage
      linkage_class[members] <- "host"
    } else {
      linkage_label[members] <- "Borg or ECE-linked"
      linkage_class[members] <- "borg_ece"
    }
  }

  # ---- Restrict to motifs well supported in this (sub)set of genomes ----
  # Keep a motif if EITHER (a) its max modification fraction and its max modified-site count
  # across the included genomes both clear the thresholds, OR (b) its max modified-site count
  # is very high (motif_min_sites_strong) regardless of fraction. Modified-site counts come
  # from the cross-profile files (motif_modified_num) so cross-detected motifs are covered.
  # Used for single-sample panels so only robust motifs are shown/labeled.
  if (motif_min_fraction > 0 || motif_min_sites > 0 || is.finite(motif_min_sites_strong)) {
    fmax <- apply(pivot_matrix, 2, max)
    if (!is.null(site_count_file) && file.exists(site_count_file)) {
      sc <- read.csv(site_count_file, stringsAsFactors = FALSE)
      sc <- sc[sc$contig %in% contig_names, ]
      site_max <- tapply(sc$nDetected, sc$motif, max, na.rm = TRUE)
      smax <- site_max[colnames(pivot_matrix)]
      smax[is.na(smax)] <- 0
    } else {
      smax <- rep(0, ncol(pivot_matrix))
    }
    crit1 <- (fmax >= motif_min_fraction) & (smax >= motif_min_sites)
    crit2 <- smax >= motif_min_sites_strong
    keep_motifs <- crit1 | crit2
    pivot_matrix <- pivot_matrix[, keep_motifs, drop = FALSE]
    cat(sprintf(paste0("Kept %d motifs of %d [frac>=%.2f & sites>=%d: %d; ",
                       "or sites>=%s: %d].\n"),
                sum(keep_motifs), length(keep_motifs), motif_min_fraction, motif_min_sites,
                sum(crit1), format(motif_min_sites_strong), sum(crit2)))
  }

  # ---- Remove redundant (containment/subset) motifs, keeping the one with MORE sites ----
  # Collapses IUPAC subset/superset and reverse-complement redundancy (aligned at the
  # modified base) to the best-supported representative (most modified sites).
  if (dedup_motifs && ncol(pivot_matrix) > 1) {
    # site weight per current motif (max modified sites across the included genomes)
    site_w <- setNames(rep(0, ncol(pivot_matrix)), colnames(pivot_matrix))
    if (!is.null(site_count_file) && file.exists(site_count_file)) {
      scd <- read.csv(site_count_file, stringsAsFactors = FALSE)
      scd <- scd[scd$contig %in% contig_names, ]
      sm <- tapply(scd$nDetected, scd$motif, max, na.rm = TRUE)
      hit <- intersect(colnames(pivot_matrix), names(sm))
      site_w[hit] <- sm[hit]
    }
    keep_ids <- dedup_motifs_by_support(colnames(pivot_matrix), site_w)
    cat(sprintf("Motif dedup (keep most-supported): %d -> %d motifs (removed %d redundant).\n",
                ncol(pivot_matrix), length(keep_ids), ncol(pivot_matrix) - length(keep_ids)))
    pivot_matrix <- pivot_matrix[, colnames(pivot_matrix) %in% keep_ids, drop = FALSE]
  }

  # ---- Drop zero-variance columns/rows (same as plot_borg_heatmap.R) ----
  col_vars <- apply(pivot_matrix, 2, var)
  valid_cols <- col_vars > 0 & !is.na(col_vars)
  if (sum(valid_cols) < 2) stop("Not enough motifs with variance for clustering.")
  pivot_matrix <- pivot_matrix[, valid_cols]
  cat(sprintf("Removed %d zero-variance motifs. Keeping %d.\n",
              sum(!valid_cols), sum(valid_cols)))

  row_vars <- apply(pivot_matrix, 1, var)
  valid_rows <- row_vars > 0 & !is.na(row_vars)
  if (sum(valid_rows) < 2) stop("Not enough contigs with variance for clustering.")
  if (sum(!valid_rows) > 0) {
    cat(sprintf("Removed %d zero-variance contigs. Keeping %d.\n",
                sum(!valid_rows), sum(valid_rows)))
    pivot_matrix  <- pivot_matrix[valid_rows, ]
    contig_names  <- contig_names[valid_rows]
    y_labels_all  <- y_labels_all[valid_rows]
    type_names    <- type_names[valid_rows]
    comp_id       <- comp_id[valid_rows]
    linkage_label <- linkage_label[valid_rows]
    linkage_class <- linkage_class[valid_rows]
    sample_names  <- sample_names[valid_rows]
  }

  # ---- Column clustering (rows are grouped by linkage group, see row_split below) ----
  col_hclust <- hclust(as.dist(1 - cor(pivot_matrix, method = "pearson")),
                       method = "ward.D2")

  # ---- Order the linkage-group slices: host-linked first (by lineage), then Borg/ECE-only
  #      linked groups, then unlinked singletons; rows within a slice cluster by profile. ----
  comp_tbl <- data.frame(comp = comp_id, cls = linkage_class, lab = linkage_label,
                         stringsAsFactors = FALSE) %>%
    group_by(comp) %>%
    summarise(cls = first(cls), lab = first(lab), size = n(), .groups = "drop")
  class_rank <- c(host = 1, borg_ece = 2, unlinked = 3)
  comp_tbl <- comp_tbl %>%
    mutate(rank = class_rank[cls]) %>%
    arrange(rank, lab, desc(size), comp)
  ordered_comps <- as.character(comp_tbl$comp)
  row_split_fac <- factor(as.character(comp_id), levels = ordered_comps)

  # ---- Colors ----
  # Type colors: reuse the paper's network palette so the figure stays consistent.
  ncdf <- read.csv(network_colors_file, stringsAsFactors = FALSE)
  type_color_map <- ncdf %>%
    filter(!is.na(genome), !is.na(Color)) %>%
    group_by(genome, Color) %>% summarise(n = n(), .groups = 'drop') %>%
    group_by(genome) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup()
  type_colors <- setNames(type_color_map$Color, type_color_map$genome)
  # Ensure every observed type has a color
  for (tp in unique(type_names)) if (is.na(type_colors[tp])) type_colors[tp] <- "#999999"
  type_colors <- type_colors[names(type_colors) %in% unique(type_names)]

  # Sample colors: reuse the Dark2 + Set1 scheme (yellow removed) from plot_borg_heatmap.R.
  unique_samples <- unique(sample_names)
  n_samples <- length(unique_samples)
  if (n_samples <= 8) {
    sample_pal <- brewer.pal(min(max(3, n_samples), 8), "Dark2")[1:n_samples]
  } else {
    base_colors <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))
    base_colors <- base_colors[base_colors != "#FFFF33"]
    sample_pal <- if (n_samples <= length(base_colors)) base_colors[1:n_samples]
                  else colorRampPalette(base_colors)(n_samples)
  }
  sample_colors <- setNames(sample_pal, unique_samples)

  # Heatmap value colormap (same yellow -> blue ramp)
  col_fun <- colorRamp2(seq(0, 1, length.out = 100),
                        colorRampPalette(c("#FFFFD4", "#C7E9B4", "#41B6C4",
                                           "#225EA8", "#081D58"))(100))

  # ---- Left annotation. The linkage grouping is conveyed by the row-split blocks and gaps,
  #      so a separate linkage-group bar is not needed. The Sample bar is dropped when the
  #      plot covers a single sample (nothing to distinguish). ----
  if (length(unique(sample_names)) > 1) {
    left_anno <- rowAnnotation(
      `Type` = type_names,
      `Sample` = sample_names,
      col = list(`Type` = type_colors, `Sample` = sample_colors),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_legend_param = list(
        `Type` = list(title = "Type", title_gp = gpar(fontsize = 9, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8)),
        `Sample` = list(title = "Sample", title_gp = gpar(fontsize = 9, fontface = "bold"),
                        labels_gp = gpar(fontsize = 7))
      ),
      width = unit(2 * 4, "mm")
    )
  } else {
    left_anno <- rowAnnotation(
      `Type` = type_names,
      col = list(`Type` = type_colors),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_legend_param = list(
        `Type` = list(title = "Type", title_gp = gpar(fontsize = 9, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8))
      ),
      width = unit(4, "mm")
    )
  }

  hm_legend <- list(title = "Modification fraction",
                    title_gp = gpar(fontsize = 9, fontface = "bold"),
                    labels_gp = gpar(fontsize = 8))

  if (!transpose) {
    # Genomes as rows (split into linkage-group blocks, separated by gaps); motifs as columns.
    ht <- Heatmap(pivot_matrix,
                  name = "Modification fraction",
                  col = col_fun,
                  row_split = row_split_fac,
                  cluster_rows = TRUE,
                  clustering_distance_rows = "pearson",
                  clustering_method_rows = "ward.D2",
                  cluster_row_slices = FALSE,
                  show_row_dend = FALSE,
                  row_gap = unit(1.2, "mm"),
                  row_title = NULL,
                  cluster_columns = col_hclust,
                  show_row_names = TRUE,
                  show_column_names = show_motif_names,
                  show_heatmap_legend = TRUE,
                  left_annotation = left_anno,
                  row_labels = y_labels_all,
                  row_names_gp = gpar(fontsize = 6),
                  row_names_max_width = unit(20, "cm"),
                  column_names_gp = gpar(fontsize = 5),
                  column_names_rot = 90,
                  column_dend_height = unit(1, "cm"),
                  heatmap_legend_param = hm_legend)
  } else {
    # Rotated layout for single-sample panels: motifs as rows (readable horizontal labels),
    # genomes as columns split into the same linkage-group blocks with gaps.
    mat_t <- t(pivot_matrix)
    # auto-size: motifs on rows (height), genomes on columns (width)
    plot_height <- max(5, nrow(mat_t) * 0.13 + 4)
    plot_width  <- max(4.5, ncol(mat_t) * 0.28 + 4)
    col_labels_short <- mapply(short_genome_label, contig_names,
                               type_names, contig_to_borg[contig_names])
    top_anno <- HeatmapAnnotation(
      `Type` = type_names,
      col = list(`Type` = type_colors),
      annotation_name_gp = gpar(fontsize = 8),
      annotation_legend_param = list(
        `Type` = list(title = "Type", title_gp = gpar(fontsize = 9, fontface = "bold"),
                      labels_gp = gpar(fontsize = 8))),
      height = unit(4, "mm"), which = "column")
    ht <- Heatmap(mat_t,
                  name = "Modification fraction",
                  col = col_fun,
                  column_split = row_split_fac,
                  cluster_columns = TRUE,
                  clustering_distance_columns = "pearson",
                  clustering_method_columns = "ward.D2",
                  cluster_column_slices = FALSE,
                  show_column_dend = FALSE,
                  column_gap = unit(1.2, "mm"),
                  column_title = NULL,
                  cluster_rows = col_hclust,
                  show_row_dend = TRUE,
                  row_dend_width = unit(1, "cm"),
                  show_row_names = show_motif_names,
                  show_column_names = TRUE,
                  show_heatmap_legend = TRUE,
                  top_annotation = top_anno,
                  column_labels = col_labels_short,
                  row_names_gp = gpar(fontsize = 6),
                  column_names_gp = gpar(fontsize = 7),
                  column_names_rot = 90,
                  heatmap_legend_param = hm_legend)
  }

  pdf(plot_name, width = plot_width, height = plot_height)
  draw(ht, merge_legend = TRUE,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
  cat(paste0("Saved annotated heatmap to: ", plot_name, "\n"))

  # ---- Source data (Nature "Source Data" convention) ----
  sourcedata <- data.frame(
    contig = contig_names,
    type = unname(type_names),
    linkage_component = unname(comp_id),
    linkage_class = unname(linkage_class),
    linkage_group = unname(linkage_label),
    slice_order = match(as.character(comp_id), ordered_comps),
    sample = unname(sample_names),
    stringsAsFactors = FALSE
  ) %>% arrange(slice_order, contig)
  sd_file <- sub("\\.pdf$", ".sourcedata.csv", plot_name)
  write.csv(sourcedata, sd_file, row.names = FALSE)
  cat(paste0("Saved source data to: ", sd_file, "\n"))

  # ---- Console sanity check on network-derived linkage grouping ----
  cat("\nNetwork linkage grouping by type (rows split into these blocks):\n")
  tab <- table(type_names, linkage_class)
  print(tab)
  cat(sprintf("\n%d linkage-group blocks total (%d multi-member groups, %d unlinked singletons).\n",
              length(ordered_comps),
              sum(comp_tbl$size >= 2), sum(comp_tbl$size < 2)))

  invisible(NULL)
}

# ---------------------------------------------------------------------------
seq_dir <- "/home/shuaiw/borg/paper/borg_data/profile4"   # input data location
out_dir <- "/home/shuaiw/MODIFI/tmp/rev_figs/borg"        # figure + source data output
cluster <- "profile"
script_dir <- tryCatch(dirname(sub("^--file=", "",
                    grep("^--file=", commandArgs(FALSE), value = TRUE))),
                    error = function(e) ".")
if (length(script_dir) == 0 || script_dir == "") script_dir <- "."

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

profile_df <- read.csv(file.path(seq_dir, paste0(cluster, "_profile_df_filtered.csv")))
profile_df <- profile_df %>% filter(Genome != "Non-Mp")

# sample column is needed by the network step in profile_borg.py but not by us here;
# we take sample from the sample.name file instead.

edge_file    <- file.path(seq_dir, "borg_motif_profile_all_network_edges.csv")
colors_file  <- file.path(script_dir, "network_colors.csv")
sample_file  <- file.path(seq_dir, paste0(cluster, "_profile_df_filtered.sample.name.csv"))
sites_file   <- file.path(seq_dir, "motif_site_counts.csv")  # built by build_motif_site_counts.py

# ---- Full figure (all samples) ----
personal_plot_annotated(
  profile_df,
  plot_name           = file.path(out_dir, "borg_profile_heatmap_annotated.pdf"),
  network_edge_file   = edge_file,
  network_colors_file = colors_file,
  sample_name_file    = sample_file,
  plot_height         = 11
)

# ---- Single-sample panel(s): one plot per sample ----
sample_map <- read.csv(sample_file, stringsAsFactors = FALSE)
all_samples <- sort(unique(profile_df %>%
  inner_join(sample_map, by = "contig") %>% pull(sample)))
for (target_sample in all_samples) {
  target_contigs <- sample_map$contig[sample_map$sample == target_sample]
  sub_df <- profile_df %>% filter(contig %in% target_contigs)
  n_ctg <- length(unique(sub_df$contig))
  if (n_ctg < 2) {
    cat(sprintf("Skipping %s: only %d contig(s).\n", target_sample, n_ctg))
    next
  }
  cat(sprintf("\n=== Single-sample panel: %s (%d contigs) ===\n", target_sample, n_ctg))
  safe <- gsub("[^A-Za-z0-9._-]", "_", target_sample)
  personal_plot_annotated(
    sub_df,
    plot_name           = file.path(out_dir, paste0("borg_profile_heatmap_annotated_", safe, ".pdf")),
    network_edge_file   = edge_file,
    network_colors_file = colors_file,
    sample_name_file    = sample_file,
    plot_height         = max(4.5, n_ctg * 0.28 + 3.5),
    motif_min_fraction     = 0.5,   # criterion 1: max fraction >= 0.5 ...
    motif_min_sites        = 100,   #             AND max modified sites >= 100
    motif_min_sites_strong = 1000,  # criterion 2: OR max modified sites >= 1000 (any fraction)
    site_count_file        = sites_file,
    show_motif_names    = TRUE,  # label motif strings (now horizontal row labels)
    dedup_motifs        = TRUE,  # collapse containment-redundant motifs (keep specific)
    transpose           = TRUE   # rotate: motifs as rows, genomes as columns
  )
}
