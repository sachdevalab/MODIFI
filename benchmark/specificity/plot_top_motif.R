library(tidyverse)
library(readr)
library(viridis)

# Set figure directory
fig_dir <- "../../tmp/figures/multi_env_linkage"

# Function to get detail taxa name from lineage
get_detail_taxa_name <- function(lineage) {
  if (lineage == "Unknown") {
    return("Unknown")
  }
  
  # Split lineage by semicolon
  levels <- str_split(lineage, ";")[[1]]
  
  # Loop through levels in reverse order
  for (level in rev(levels)) {
    level_trimmed <- str_trim(level)
    if (nchar(level_trimmed) > 3) {
      return(level_trimmed)
    }
  }
  
  # If no level with > 3 characters found, return first level
  return(str_trim(levels[1]))
}

# Read data
df_all_data <- read_csv(file.path(fig_dir, "motif_num_all_samples.csv"))

# Extract species from lineage using the get_detail_taxa_name function
# Remove p__ prefix from phylum names
df_all_data <- df_all_data %>%
  mutate(
    species = sapply(lineage, get_detail_taxa_name),
    phylum = str_replace(phylum, "^p__", "")
  )

# Get top 20 genomes by motif number
top20_df <- df_all_data %>%
  arrange(desc(motif_num)) %>%
  head(20)

# Define phylum color palette
PHYLUM_COLORS <- c(
  "Pseudomonadota" = "#d8b365",
  "Bacillota" = "#f46d43",
  "Bacillota_A" = "#8da0cb",
  "Desulfobacterota" = "#fc8d62",
  "Actinomycetota" = "#66c2a5",
  "Bacteroidota" = "#e78ac3",
  "Campylobacterota" = "#a6d854",
  "Acidobacteriota" = "#1b9e77",
  "Verrucomicrobiota" = "#7570b3",
  "Chloroflexota" = "#e7298a",
  # Archaea phyla
  "Thermoproteota" = "#b15928",
  "Thermoplasmatota" = "#fdbf6f",
  "Halobacteriota" = "#ff7f00",
  "Aenigmatarchaeota" = "#cab2d6",
  "Nanoarchaeota" = "#6a3d9a",
  "Hadarchaeota" = "#fb9a99",
  "Asgardarchaeota" = "#e31a1c",
  "Micrarchaeota" = "#a6cee3",
  "Others" = "#e6e6e6"
)

# Define environment color palette
ENV_COLORS <- c(
  "mock" = "#e41a1c",
  "mice_gut" = "#377eb8",
  "sugarcane" = "#4daf4a",
  "infant_gut" = "#984ea3",
  "cow_rumen" = "#ff7f00",
  "cow_bioreactor" = "#ffff33",
  "adult_gut" = "#a65628",
  "ocean" = "#f781bf",
  "soil" = "#999999"
)

# Create combined labels that include species and contig (without n= prefix)
# Map phylum to color category (use "Others" for phyla not in the palette)
top20_df <- top20_df %>%
  mutate(
    combined_label = paste0(species, "\n(", contig, ")"),
    phylum_color = ifelse(phylum %in% names(PHYLUM_COLORS), phylum, "Others")
  )

# Create horizontal barplot with motif numbers displayed in the bars
p <- ggplot(top20_df, aes(x = motif_num, y = reorder(combined_label, motif_num))) +
  geom_col(aes(fill = phylum_color), show.legend = FALSE) +
  # Add environment color indicator on the right side of bars
  geom_segment(aes(x = motif_num + 0.5, xend = motif_num + 0.5, 
                   y = as.numeric(reorder(combined_label, motif_num)) - 0.4,
                   yend = as.numeric(reorder(combined_label, motif_num)) + 0.4,
                   color = environment), size = 3) +
  geom_text(aes(label = motif_num), hjust = 1.2, color = "white", size = 4, fontface = "bold") +
  scale_fill_manual(values = PHYLUM_COLORS) +
  scale_color_manual(values = ENV_COLORS, name = "Environment") +
  labs(x = "Number of Motifs",
       y = "") +
    #    y = "Taxa (contig)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        legend.position = "top",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10))

# Save plot as PDF
ggsave(file.path(fig_dir, "top20_genome_motif.pdf"), 
       plot = p, 
       width = 6, 
       height = 10)