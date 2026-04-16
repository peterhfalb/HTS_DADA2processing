# 06_combine_otu_taxonomy.R
# Combines mumu-curated OTU abundance table with taxonomy/bootstrap assignments.
# Renames sample columns from bioinformatics IDs (S001...) back to original names.
#
# Usage: Rscript 06_combine_otu_taxonomy.R <proj_name> <primer_set> <otu_table> <taxonomy_combined> <map_file> <output_file>

library(dplyr)
library(readr)
library(tibble)

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript 06_combine_otu_taxonomy.R <proj_name> <primer_set> <otu_table> <taxonomy_combined> <map_file> <output_file>")
}

PROJ              <- args[1]
PRIMER_SET        <- args[2]
otu_table_file    <- args[3]
tax_combined_file <- args[4]
map_file          <- args[5]
output_file       <- args[6]

cat("Project name:  ", PROJ, "\n")
cat("Primer set:    ", PRIMER_SET, "\n")
cat("OTU table:     ", otu_table_file, "\n")
cat("Taxonomy file: ", tax_combined_file, "\n")
cat("Map file:      ", map_file, "\n")
cat("Output file:   ", output_file, "\n")

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------

cat("\nLoading taxonomy table...\n")
tax <- read_tsv(tax_combined_file, show_col_types = FALSE)
cat("  OTUs in taxonomy table: ", nrow(tax), "\n")
cat("  Taxonomy columns: ", paste(names(tax), collapse = ", "), "\n")

cat("\nLoading OTU abundance table...\n")
otu <- read_tsv(otu_table_file, show_col_types = FALSE)
cat("  OTUs in mumu-curated table: ", nrow(otu), "\n")

cat("\nLoading sample map file...\n")
map <- read_csv(map_file, show_col_types = FALSE)
cat("  Samples in map: ", nrow(map), "\n")
print(map)

# ------------------------------------------------------------------------------
# Rename OTU table sample columns from bioinformatics IDs to original names
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Rename OTU table sample columns from bioinformatics IDs to original names
# (only if columns are still in S001/S002 format)
# ------------------------------------------------------------------------------

otu_id_col  <- names(otu)[1]
sample_cols <- names(otu)[-1]

# Check if renaming is needed (columns look like S001, S002, ..., S10000+)
needs_rename <- all(grepl("^S\\d+$", sample_cols))

if (needs_rename) {
  rename_map <- setNames(map$Sample_name, map$BioInf_Sample)
  missing <- setdiff(sample_cols, names(rename_map))
  if (length(missing) > 0) {
    warning("Some sample columns not found in Map_file.csv: ", paste(missing, collapse = ", "))
  }
  names(otu)[-1] <- rename_map[sample_cols]
  cat("\nSample columns renamed to original names.\n")
} else {
  cat("\nSample columns already have original names, skipping rename.\n")
}

# ------------------------------------------------------------------------------
# Join OTU table with taxonomy
# # OTU ID column in mumu table may have ;size= suffixes stripped already,
# but taxonomy table uses the clean centroid ID — match on first column
# ------------------------------------------------------------------------------

# Ensure OTU ID column name is consistent for joining
names(otu)[1]  <- "OTU_ID"
names(tax)[1]  <- "OTU_ID"

combined <- otu %>%
  inner_join(tax, by = "OTU_ID")

cat("\nOTUs after joining OTU table with taxonomy: ", nrow(combined), "\n")

if (nrow(combined) < nrow(otu)) {
  cat("  WARNING: ", nrow(otu) - nrow(combined),
      " OTUs from abundance table had no taxonomy match and were dropped.\n")
}

# ------------------------------------------------------------------------------
# Write output
# Column order: OTU_ID | sample abundance columns | taxonomy + bootstrap columns
# ------------------------------------------------------------------------------

sample_col_names <- map$Sample_name
tax_col_names    <- setdiff(names(combined), c("OTU_ID", sample_col_names))

final <- combined %>%
  select(OTU_ID, all_of(sample_col_names), all_of(tax_col_names))

write_tsv(final, output_file)

cat("\nFinal output written to: ", output_file, "\n")
cat("  OTUs:    ", nrow(final), "\n")
cat("  Samples: ", length(sample_col_names), "\n")
cat("  Columns: OTU_ID +", length(sample_col_names), "sample columns +",
    length(tax_col_names), "taxonomy/bootstrap columns\n")