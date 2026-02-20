fig_dir <- "../../tmp/figures/multi_env_linkage/"
df_all_data <- read.csv(paste0(fig_dir, "motif_num_all_samples.csv"))

# Count genomes with 0 motifs per phylum
zero_motif_genomes <- df_all_data[df_all_data$motif_num == 0, ]
zero_counts <- table(zero_motif_genomes$phylum)
zero_counts_sorted <- sort(zero_counts, decreasing = TRUE)

cat("\nNumber of genomes with 0 motifs per phylum:\n")
print(zero_counts_sorted)
cat("\nTotal genomes with 0 motifs:", nrow(zero_motif_genomes), "\n")

# Extract family from lineage and count genomes with 0 motifs per family
extract_family <- function(lineage) {
  parts <- strsplit(as.character(lineage), ";")[[1]]
  family_part <- parts[grepl("^f__", parts)]
  if (length(family_part) > 0) {
    return(family_part[1])
  } else {
    return(NA)
  }
}

zero_motif_genomes$family <- sapply(zero_motif_genomes$lineage, extract_family)
zero_counts_family <- table(zero_motif_genomes$family)
zero_counts_family_sorted <- sort(zero_counts_family, decreasing = TRUE)

cat("\n\nNumber of genomes with 0 motifs per family:\n")
print(zero_counts_family_sorted)

# Output RM_num for each genome with 0 motifs
cat("\n\nRM_num for each genome with 0 motifs:\n")
zero_motif_output <- zero_motif_genomes[, c("contig", "phylum", "family", "RM_num")]
print(zero_motif_output, row.names = FALSE)

# Summary statistics of RM_num for genomes with 0 motifs
cat("\n\nRM_num summary for genomes with 0 motifs:\n")
cat("Mean:", mean(zero_motif_genomes$RM_num), "\n")
cat("Median:", median(zero_motif_genomes$RM_num), "\n")
cat("Min:", min(zero_motif_genomes$RM_num), "\n")
cat("Max:", max(zero_motif_genomes$RM_num), "\n")
cat("Genomes with RM_num = 0:", sum(zero_motif_genomes$RM_num == 0), "\n")
cat("Genomes with RM_num > 0:", sum(zero_motif_genomes$RM_num > 0), "\n")