# 03_build_otu_table.R
# Builds an OTU abundance table by mapping ASVs to OTU centroids via the VSEARCH .uc file
# then summing per-sample abundances from the original ASV abundance table
# Usage: Rscript 03_build_otu_table.R <uc_file> <asv_abundance_table> <output_file>

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(purrr)
library(stringr)
# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 03_build_otu_table.R <uc_file> <asv_abundance_table> <output_file>")
}

UC_FILE    <- args[1]
ASV_FILE   <- args[2]
OUTPUT     <- args[3]

cat("UC file:           ", UC_FILE, "\n")
cat("ASV abundance file:", ASV_FILE, "\n")
cat("Output file:       ", OUTPUT, "\n")

# ------------------------------------------------------------------------------
# Parse .uc file
# UC format columns:
#   1=type, 2=cluster_nr, 3=size, 4=pct_id, 5=strand,
#   6=query_start, 7=seed_start, 8=alignment, 9=query_label, 10=target_label
# We want H (hit) records: col 9 = ASV (query), col 10 = OTU centroid (target)
# S (seed) records are the centroids themselves — they map to themselves
# ------------------------------------------------------------------------------

cat("Parsing .uc file...\n")

uc_raw <- read_tsv(UC_FILE,
                   col_names = c("type", "cluster_nr", "size", "pct_id", "strand",
                                 "query_start", "seed_start", "alignment",
                                 "query_label", "target_label"),
                   col_types = cols(.default = "c"),
                   comment = "#")

# H records: ASV matched to a centroid
hits <- uc_raw %>%
  filter(type == "H") %>%
  select(query_label, target_label) %>%
  rename(SeqID = query_label, OTUId = target_label)

# S records: centroids map to themselves
seeds <- uc_raw %>%
  filter(type == "S") %>%
  select(query_label) %>%
  mutate(OTUId = query_label) %>%
  rename(SeqID = query_label)

# Combine: every ASV now has an OTU assignment
asv_to_otu <- bind_rows(hits, seeds)

# Strip size annotations from labels if present (e.g. SeqXXXXXX;size=100 -> SeqXXXXXX)
asv_to_otu <- asv_to_otu %>%
  mutate(
    SeqID = sub(";.*", "", SeqID),
    OTUId = sub(";.*", "", OTUId)
  )

cat("ASVs mapped to OTUs:", nrow(asv_to_otu), "\n")
cat("Unique OTUs:        ", n_distinct(asv_to_otu$OTUId), "\n")

# ------------------------------------------------------------------------------
# Load ASV abundance table
# ------------------------------------------------------------------------------

cat("Loading ASV abundance table...\n")

asv_abund <- read_tsv(ASV_FILE, col_types = cols(.default = "d", SeqID = "c"))

cat("ASVs in abundance table:", nrow(asv_abund), "\n")

# ------------------------------------------------------------------------------
# Build OTU table by summing ASV abundances per OTU
# ------------------------------------------------------------------------------

cat("Building OTU table...\n")

otu_table <- asv_abund %>%
  left_join(asv_to_otu, by = "SeqID") %>%
  filter(!is.na(OTUId)) %>%       # drop ASVs that didn't map (filtered out during clustering)
  group_by(OTUId) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  arrange(desc(rowSums(across(where(is.numeric)))))

cat("OTUs in final table:", nrow(otu_table), "\n")

# ------------------------------------------------------------------------------
# Write output
# ------------------------------------------------------------------------------

write.table(otu_table, OUTPUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat("OTU table written to:", OUTPUT, "\n")