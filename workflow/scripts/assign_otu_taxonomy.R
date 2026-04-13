# RDP Taxonomy Assignment via DADA2
# Usage: Rscript 05_assign_taxonomy_rdp.R <fasta> <database> <output> <threads>

options(repos = c(CRAN = "https://cloud.r-project.org"))  # set CRAN mirror

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(purrr)
library(stringr)
library(dada2)
library(Biostrings)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 05_assign_taxonomy_rdp.R <fasta_path> <db_path> <output_path> <threads>")
}

FASTA_PATH   <- args[1]
DB_PATH      <- args[2]
OUTPUT_PATH  <- args[3]
PRIMER_SET   <- args[4]
DB_NAME    <- ifelse(length(args) >= 5, args[5], "")

cat("Input FASTA:  ", FASTA_PATH, "\n")
cat("Database:     ", DB_PATH, "\n")
cat("Output:       ", OUTPUT_PATH, "\n")

# ------------------------------------------------------------------------------
# Read sequences from FASTA
# ------------------------------------------------------------------------------

seqs <- readDNAStringSet(FASTA_PATH)
seq_ids <- names(seqs)
sequences <- as.character(seqs)

# Strip size annotations (e.g., Seq000001;size=100 -> Seq000001)
seq_ids <- sub(";.*", "", seq_ids)

cat("Sequences to classify:", length(sequences), "\n")

# ------------------------------------------------------------------------------
# Assign taxonomy with bootstrap values
# ------------------------------------------------------------------------------

cat("Running assignTaxonomy...\n")

# Set taxonomy levels based on primer set and database
# These must exactly match the number of semicolon-delimited levels in the database headers
if (PRIMER_SET == "18S-V4" && DB_NAME == "PR2") {
  # PR2 uses a non-standard 9-level hierarchy
  tax_levels <- c("Domain", "Supergroup", "Division", "Subdivision",
                  "Class", "Order", "Family", "Genus", "Species")
} else {
  # UNITE, SILVA, EukaryomeITS, EukaryomeSSU, MaarjAM — standard 7-level hierarchy
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
}

tax <- assignTaxonomy(
  sequences,
  DB_PATH,
  taxLevels        = tax_levels,
  multithread      = TRUE,
  outputBootstraps = TRUE,
  tryRC            = grepl("^18S", PRIMER_SET),  # Use reverse complement for 18S primers
  verbose          = TRUE
)

# ------------------------------------------------------------------------------
# Build tables
# ------------------------------------------------------------------------------

tax_table  <- data.frame(OTUId = seq_ids, Sequence = sequences, tax$tax,  stringsAsFactors = FALSE)
boot_table <- data.frame(OTUId = seq_ids, Sequence = sequences, tax$boot, stringsAsFactors = FALSE)
colnames(boot_table)[-(1:2)] <- paste0(colnames(boot_table)[-(1:2)], "_bootstrap")

# Combined table: OTUId, Sequence, taxonomy columns, then bootstrap columns
combined_table <- cbind(
  tax_table,
  boot_table[, -(1:2)]  # drop OTUId and Sequence since already present
)

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

write.table(tax_table,
            OUTPUT_PATH,
            sep = "\t", quote = FALSE, row.names = FALSE)

boot_path <- sub("\\.txt$", "_bootstraps.txt", OUTPUT_PATH)
write.table(boot_table,
            boot_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

combined_path <- sub("\\.txt$", "_combined.txt", OUTPUT_PATH)
write.table(combined_table,
            combined_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Taxonomy assigned and written to:        ", OUTPUT_PATH, "\n")
cat("Bootstrap values written to:             ", boot_path, "\n")
cat("Combined taxonomy+bootstraps written to: ", combined_path, "\n")