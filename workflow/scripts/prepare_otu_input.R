#!/usr/bin/env Rscript
# prepare_otu_input.R
# Prepares a single combined FASTA and abundance table from a combined ASV+taxonomy table
# Usage: Rscript prepare_otu_input.R <asv_table_path> <output_dir>

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(purrr)
library(stringr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript prepare_otu_input.R <asv_table_path> <output_dir>")
}

ASVTABLE_PATH <- args[1]
WORK_DIR      <- args[2]

cat("ASV table path:", ASVTABLE_PATH, "\n")
cat("Working dir:   ", WORK_DIR, "\n")

setwd(WORK_DIR)

# Load and parse ASV table
fullTable <- read_tsv(ASVTABLE_PATH, show_col_types = FALSE)

# Remove taxonomy and bootstrap columns by name
# Covers all primer sets: ITS, 16S, 18S-V4, 18S-AMF
tax_ranks <- c(
  "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species",
  "Domain", "Supergroup", "Division", "Subdivision",
  "NA"
)

# Build pattern that matches exact rank names AND deduplicated variants (e.g. Kingdom...196)
tax_pattern <- paste0("^(", paste(tax_ranks, collapse = "|"), ")(\\.\\.\\.[0-9]+)?$")

keep <- !grepl(tax_pattern, names(fullTable), ignore.case = FALSE)
keep[1] <- TRUE  # always keep ASV column (sequences)

asvT <- fullTable[, keep]
names(asvT)[1] <- "ASV"

cat("Columns retained (ASV + sample counts): ", ncol(asvT), "\n")

n_samples <- ncol(asvT)  # includes ASV column

# Add total reads and sort descending
asvT <- asvT %>%
  mutate(Total = rowSums(across(-1))) %>%  # all columns except the first (ASV)
  arrange(desc(Total))

# Filter out ASVs with zero total reads
asvT <- asvT %>% filter(Total > 0)
cat("ASVs after removing zero-abundance entries:", nrow(asvT), "\n")

cat("ASVs loaded:  ", nrow(asvT), "\n")
cat("Samples:      ", n_samples - 1, "\n")

# Extract sequences
sequences <- asvT$ASV

# Assign SeqXXXXXX IDs
seq_ids <- paste0("Seq", str_pad(seq_len(nrow(asvT)), width = 6, side = "left", pad = "0"))

# Build and save sample mapping file
sample_names <- names(asvT)[2:n_samples]

Map_file <- data.frame(
  BioInf_Sample = paste0("S", str_pad(seq_along(sample_names), width = 3, side = "left", pad = "0")),
  Sample_name   = sample_names,
  stringsAsFactors = FALSE
)

cat("Sample mapping:\n")
print(Map_file)

write.csv(Map_file, "Map_file.csv", row.names = FALSE, quote = FALSE)
cat("Map file written to: Map_file.csv\n")

# Save abundance table with SeqXXXXXX row IDs and original sample column names
# This is used later by build_otu_table.R to reconstruct OTU abundances
abund_table <- asvT[, 2:n_samples]
rownames(abund_table) <- seq_ids

write.table(
  cbind(SeqID = seq_ids, abund_table),
  "asv_abundance_table.txt",
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("ASV abundance table written to: asv_abundance_table.txt\n")

# Write single combined FASTA for pipeline input
# Headers: >SeqXXXXXX;size=<total_reads>
# size= is required by VSEARCH for abundance-aware dereplication and clustering
fasta_lines <- as.vector(rbind(
  paste0(">", seq_ids, ";size=", asvT$Total),
  sequences
))

writeLines(fasta_lines, "Centroid.fasta")
cat("Combined FASTA written to: Centroid.fasta (", length(sequences), " sequences)\n")
