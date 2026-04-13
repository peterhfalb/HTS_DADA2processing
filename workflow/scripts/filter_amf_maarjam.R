# AMF Dataset Filtering: MaarjAM vs EukaryomeSSU
#
# Purpose: For 18S-AMF datasets, this script performs a secondary taxonomy assignment
# using EukaryomeSSU (broad eukaryote 18S reference) and filters the MaarjAM assignments
# to retain only those sequences detected as Mucoromycota by EukaryomeSSU (with bootstrap >= 50).
#
# This dual-assignment approach:
# - Uses MaarjAM for detailed AMF identification (species-level resolution)
# - Uses EukaryomeSSU as a quality control to exclude non-Fungi and non-Mucoromycota sequences
#
# Output: Two final OTU+taxonomy files for user to choose from:
#   1. Unfiltered (standard MaarjAM)
#   2. Filtered (Mucoromycota-only, SILVA-validated)
#
# Usage: Rscript 07_filter_amf_maarjam.R <project_name> <maarjam_combined_tax> \
#          <fasta_path> <silva_db> <otu_abundance_table> <map_file> \
#          <output_unfiltered> <output_filtered> <summary_report>

options(repos = c(CRAN = "https://cloud.r-project.org"))

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(purrr)
library(stringr)
library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript 07_filter_amf_maarjam.R <proj_name> <maarjam_combined_tax> ",
       "<fasta_path> <silva_db> <otu_abundance_table> <map_file> ",
       "<output_unfiltered> <output_filtered> <summary_report>")
}

PROJ             <- args[1]
MAARJAM_TAX_PATH <- args[2]
FASTA_PATH       <- args[3]
EUKARYOME_DB_PATH <- args[4]
OTU_TABLE_PATH   <- args[5]
MAP_FILE_PATH    <- args[6]
OUT_UNFILTERED   <- args[7]
OUT_FILTERED     <- args[8]
SUMMARY_PATH     <- args[9]

# Derive taxonomy folder from summary_path (which is in the taxonomy folder)
TAXONOMY_DIR <- dirname(SUMMARY_PATH)
EUKARYOME_TAX_OUTPUT <- file.path(TAXONOMY_DIR, "Taxonomy_EukaryomeSSU_validation_combined.txt")

cat("Project name:                    ", PROJ, "\n")
cat("MaarjAM taxonomy:                ", MAARJAM_TAX_PATH, "\n")
cat("Input FASTA:                     ", FASTA_PATH, "\n")
cat("EukaryomeSSU reference DB:       ", EUKARYOME_DB_PATH, "\n")
cat("OTU abundance table:             ", OTU_TABLE_PATH, "\n")
cat("Map file:                        ", MAP_FILE_PATH, "\n")
cat("Output (unfiltered):             ", OUT_UNFILTERED, "\n")
cat("Output (filtered):               ", OUT_FILTERED, "\n")
cat("Summary report:                  ", SUMMARY_PATH, "\n")
cat("EukaryomeSSU validation output:  ", EUKARYOME_TAX_OUTPUT, "\n\n")

# ------------------------------------------------------------------------------
# Read input files
# ------------------------------------------------------------------------------

cat("Reading input files...\n")

# Read MaarjAM taxonomy assignments (combined with bootstraps)
maarjam_tax <- read_delim(MAARJAM_TAX_PATH, delim = "\t", col_types = cols(.default = "c"))

# Read sequences
seqs <- readDNAStringSet(FASTA_PATH)
seq_ids <- names(seqs)
sequences <- as.character(seqs)

cat("Sequences to filter:    ", length(sequences), "\n")

# Read OTU abundance table
otu_table <- read_delim(OTU_TABLE_PATH, delim = "\t", col_types = cols(.default = "c"))

# Read map file (sample mapping)
map_file <- read_csv(MAP_FILE_PATH, col_types = cols(.default = "c"))

# ------------------------------------------------------------------------------
# Run SILVA taxonomy assignment on same centroid sequences
# ------------------------------------------------------------------------------

cat("\nRunning EukaryomeSSU reference taxonomy assignment...\n")

eukaryome_tax <- assignTaxonomy(
  sequences,
  EUKARYOME_DB_PATH,
  taxLevels        = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
  multithread      = TRUE,
  outputBootstraps = TRUE,
  tryRC            = TRUE,
  verbose          = TRUE
)

# Extract EukaryomeSSU assignments and bootstraps
eukaryome_table  <- data.frame(OTUId = seq_ids, eukaryome_tax$tax, stringsAsFactors = FALSE)
eukaryome_boot   <- data.frame(OTUId = seq_ids, eukaryome_tax$boot, stringsAsFactors = FALSE)
colnames(eukaryome_boot)[-(1)] <- paste0(colnames(eukaryome_boot)[-(1)], "_bootstrap")

# Build combined EukaryomeSSU table (for reference/validation)
eukaryome_combined <- cbind(
  eukaryome_table,
  eukaryome_boot[, -1]  # drop OTUId since already present
)

# Build filtering summary table
filter_summary <- data.frame(
  OTUId = seq_ids,
  Sequence = sequences,
  maarjam_kingdom = maarjam_tax$Kingdom[match(seq_ids, maarjam_tax$OTUId)],
  maarjam_phylum = maarjam_tax$Phylum[match(seq_ids, maarjam_tax$OTUId)],
  maarjam_kingdom_boot = maarjam_tax$Kingdom_bootstrap[match(seq_ids, maarjam_tax$OTUId)],
  eukaryome_kingdom = eukaryome_table$Kingdom,
  eukaryome_phylum = eukaryome_table$Phylum,
  eukaryome_phylum_boot = eukaryome_boot$Phylum_bootstrap,
  stringsAsFactors = FALSE
)

# Determine filter decision
filter_summary <- filter_summary %>%
  mutate(
    Eukaryome_Phylum_Numeric = as.numeric(eukaryome_phylum_boot),
    Filter_Decision = case_when(
      is.na(eukaryome_kingdom) | eukaryome_kingdom == "" ~ "REMOVED",
      eukaryome_kingdom != "Fungi" ~ "REMOVED",
      is.na(eukaryome_phylum) | eukaryome_phylum == "" ~ "REMOVED",
      Eukaryome_Phylum_Numeric < 50 ~ "REMOVED",
      eukaryome_phylum != "Mucoromycota" ~ "REMOVED",
      TRUE ~ "KEPT"
    ),
    Reason = case_when(
      is.na(eukaryome_kingdom) | eukaryome_kingdom == "" ~ "Unclassified_at_Kingdom",
      eukaryome_kingdom != "Fungi" ~ paste0("Non_Fungi_detected_as_", eukaryome_kingdom),
      is.na(eukaryome_phylum) | eukaryome_phylum == "" ~ "Unclassified_at_Phylum",
      Eukaryome_Phylum_Numeric < 50 ~ "Low_Phylum_bootstrap",
      eukaryome_phylum != "Mucoromycota" ~ paste0("Non_Mucoromycota_detected_as_", eukaryome_phylum),
      TRUE ~ "Mucoromycota_Fungi"
    )
  ) %>%
  select(-Eukaryome_Phylum_Numeric)

cat("\nFiltering Summary:\n")
kept_count <- sum(filter_summary$Filter_Decision == "KEPT")
removed_count <- sum(filter_summary$Filter_Decision == "REMOVED")
cat("  Total sequences processed:     ", nrow(filter_summary), "\n")
cat("  Sequences KEPT (Mucoromycota): ", kept_count, "\n")
cat("  Sequences REMOVED:             ", removed_count, "\n")

# Breakdown of removals
if (removed_count > 0) {
  cat("\nRemoval Reason Breakdown:\n")
  removal_breakdown <- filter_summary %>%
    filter(Filter_Decision == "REMOVED") %>%
    count(Reason) %>%
    arrange(desc(n))
  for (i in seq_len(nrow(removal_breakdown))) {
    cat("  ", removal_breakdown$Reason[i], ": ", removal_breakdown$n[i], "\n")
  }
}

# Bootstrap distribution for kept sequences
if (kept_count > 0) {
  cat("\nEukaryomeSSU Phylum (Mucoromycota) Bootstrap for KEPT sequences:\n")
  kept_phylum_boots <- as.numeric(filter_summary$eukaryome_phylum_boot[filter_summary$Filter_Decision == "KEPT"])
  cat("  Mean:   ", round(mean(kept_phylum_boots, na.rm = TRUE), 2), "\n")
  cat("  Median: ", round(median(kept_phylum_boots, na.rm = TRUE), 2), "\n")
  cat("  Min:    ", round(min(kept_phylum_boots, na.rm = TRUE), 2), "\n")
  cat("  Max:    ", round(max(kept_phylum_boots, na.rm = TRUE), 2), "\n")
}

# Write detailed filtering summary to TSV
cat("\nWriting detailed filtering summary...\n")
write_delim(filter_summary, SUMMARY_PATH, delim = "\t")

# Write EukaryomeSSU taxonomy for reference/validation
cat("Writing EukaryomeSSU validation taxonomy table...\n")
write.table(eukaryome_combined,
            EUKARYOME_TAX_OUTPUT,
            sep = "\t", quote = FALSE, row.names = FALSE)

# ------------------------------------------------------------------------------
# Process OTU tables: Create filtered and unfiltered versions
# ------------------------------------------------------------------------------

cat("\nProcessing OTU abundance tables...\n")

# Get OTU IDs to keep (based on filter decision)
kept_otus <- filter_summary$OTUId[filter_summary$Filter_Decision == "KEPT"]

# Create unfiltered version (standard MaarjAM, all OTUs)
otu_unfiltered <- otu_table

# Create filtered version (only kept OTUs)
otu_filtered <- otu_table %>%
  filter(OTUId %in% kept_otus)

cat("OTU counts:\n")
cat("  Unfiltered (all OTUs):        ", nrow(otu_unfiltered), "\n")
cat("  Filtered (Mucoromycota only): ", nrow(otu_filtered), "\n")
cat("  OTUs removed by filter:       ", nrow(otu_unfiltered) - nrow(otu_filtered), "\n")

# Get corresponding taxonomy rows
maarjam_unfiltered <- maarjam_tax
maarjam_filtered <- maarjam_tax %>%
  filter(OTUId %in% kept_otus)

# Combine OTU tables with MaarjAM taxonomy for both versions
# Match sample columns from original OTU table
sample_cols <- setdiff(colnames(otu_table), "OTUId")

final_unfiltered <- otu_unfiltered %>%
  left_join(maarjam_unfiltered, by = "OTUId") %>%
  select(OTUId, all_of(sample_cols), everything())

final_filtered <- otu_filtered %>%
  left_join(maarjam_filtered, by = "OTUId") %>%
  select(OTUId, all_of(sample_cols), everything())

# Write final output files
cat("\nWriting final OTU+taxonomy tables...\n")
write.table(final_unfiltered,
            OUT_UNFILTERED,
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(final_filtered,
            OUT_FILTERED,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nFiltering complete!\n")
cat("Outputs written to:\n")
cat("  Unfiltered OTU table: ", OUT_UNFILTERED, "\n")
cat("  Filtered OTU table:   ", OUT_FILTERED, "\n")
cat("  Filtering summary:    ", SUMMARY_PATH, "\n")
cat("  EukaryomeSSU validation taxonomy: ", EUKARYOME_TAX_OUTPUT, "\n")
