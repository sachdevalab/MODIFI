# filter_bam.R

# Load libraries
library(optparse)
library(data.table)
setDTthreads(threads = 48)
library(tidyverse)


# #### Define script options ####
# option_list = list(make_option(c("--prefix"), type = "character", default = NULL, 
#                                help = "File name prefix", metavar = "character")) 

# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)

# # Handle NULL arguments
# if (is.null(opt$prefix)) {
#   print_help(opt_parser)
#   stop("Supply file prefix", call. = FALSE)
# }


# #### Set file and QC specifications ####
# prefix <- opt$prefix  # Set file prefix

# metagenome <- gsub("_HiC.*", "", prefix)  # Extract metagenome

# Set working directory
# setwd(paste0("/groups/diamond/projects/animal/rumen/RuReacBro20203/HiC/", metagenome, "_HiC/", metagenome, "_HiC_nr_bins_circular_elements/"))
setwd("/home/shuaiw/borg/paper/run2/cow_bioreactor_1/hic")
prefix="/home/shuaiw/borg/paper/run2/cow_bioreactor_1/hic/cow_bioreactor_1_hic_markdup"

# Set QC thresholds
mapq_thresh <- 30  # Set MAPQ threshold
pident_thresh <- 95  # Set percent sequence identity threshold
align_length_thresh <- 50  # Set alignment length threshold



#### Import and format data ####
# Import and format BAM text file
print("Importing BAM...")
bam <- fread(paste0(prefix, "_align.tsv"), sep = "\t", fill = TRUE, header = FALSE, nrows = 1000)
setnames(bam,
         c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
         c("QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "MRNM", "MPOS", "ISIZE", "SEQ", "QUAL", "NM"))

# Warning if BAM text file contains non-alignment rows
if (any(!grepl("^A\\d", bam$QNAME))) {
  stop("Not all rows begin with read names")
}


#### Parse alignment information ####
# Replace '=' in MRNM with RNAME
bam[, MRNM := ifelse(MRNM == "=", RNAME, MRNM)]

# Calculate read length and sequence identity
print("Calculating read length and sequence identity...")
bam[, seq_length := nchar(SEQ)]
bam[, NM := as.numeric(str_remove(NM, "NM:i:"))]
bam[, pident := (seq_length - NM)/seq_length * 100]

# Calculate read alignment length in CIGAR strings
print("Calculating read alignment length in CIGAR strings...")
cigar_align_lengths <- bam %>%
  distinct(CIGAR) %>%
  mutate(cigar_part = str_extract_all(CIGAR, "\\D")) %>%
  mutate(cigar_part_length = str_split(str_remove(CIGAR, "\\D$"), "\\D")) %>%
  unnest(cols = c(cigar_part, cigar_part_length)) %>%
  mutate(cigar_part_length = as.numeric(cigar_part_length)) %>%
  filter(cigar_part == "M") %>%
  group_by(CIGAR) %>%
  summarize(align_length = sum(cigar_part_length)) %>%
  ungroup() %>%
  select(CIGAR, align_length)

bam <- bam %>%
  left_join(cigar_align_lengths, by = "CIGAR")

# Identify duplicate alignments
print("Identifying duplicates...")
bam[, dup := ifelse((FLAG - 1024) > 0, TRUE, FALSE)]


#### Summarize paired-end reads ####
# Calculate QC criteria for read pairs
print("Calculating QC criteria for read pairs...")
bam_pe <- bam[, .(nreads = .N,
                  min_mapq = min(MAPQ),
                  min_seq_length = min(seq_length),
                  min_pident = min(pident),
                  min_align_length = min(align_length)), 
              by = QNAME]


#### Filter read pairs ####
# Retain read pairs based on both reads passing QC criteria
print("Filtering read pairs...")
bam_pe_filt <- bam_pe[nreads == 2 &
                        min_mapq >= mapq_thresh &
                        min_pident >= pident_thresh &
                        min_align_length >= align_length_thresh]

bam_filt <- bam[QNAME %in% bam_pe_filt$QNAME]

# Write filtered reads (including duplicates) to file
unique(bam_filt[, .(QNAME)])[, fwrite(.SD, file = paste0(prefix, "_reads.txt"), col.names = FALSE)]

# Write filtered reads (without duplicates) to file
bam_filt %>%
  filter(!dup) %>%
  distinct(QNAME) %>%
  write_tsv(paste0(prefix, "_reads_rmdup.txt"), col_names = FALSE)


#### Calculate number of reads mapped to scaffolds ####
# Calculate number of reads mapped between scaffold pairs
print("Calculating number of mapped reads between scaffold pairs...")
scaffold_pairs_nreads <- bam_filt[, .(scaffold1 = pmin(RNAME, MRNM),
                                      scaffold2 = pmax(RNAME, MRNM))] %>%
  .[, .(nreads = .N), by = .(scaffold1, scaffold2)]

# Calculate number of reads mapped between scaffold pairs without duplicates
scaffold_pairs_nreads_rmdups <- bam_filt %>%
  filter(!dup) %>%
  mutate(scaffold1 = pmin(RNAME, MRNM),
         scaffold2 = pmax(RNAME, MRNM)) %>%
  group_by(scaffold1, scaffold2) %>%
  count(name = "nreads_rmdup") %>%
  ungroup()

# Combine number of reads mapped between scaffold pairs with and without duplicates
scaffold_pairs_nreads_all <- scaffold_pairs_nreads %>%
  full_join(scaffold_pairs_nreads_rmdups, by = c("scaffold1", "scaffold2")) %>%
  mutate(nreads_rmdup = case_when(is.na(nreads_rmdup) ~ 0,
                                  TRUE ~ nreads_rmdup))

# Write number of reads mapped between scaffold pairs to file
print("Writing mapped read counts between scaffold pairs...")
fwrite(scaffold_pairs_nreads_all, file = paste0(prefix, "_scaffold_pairs_nreads.tsv"), sep = "\t", col.names = FALSE)

print("All done!")
