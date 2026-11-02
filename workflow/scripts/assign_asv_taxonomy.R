#!/usr/bin/env Rscript
# DADA2 RDP taxonomy assignment for ASVs, and the final ASV tables
# Usage: Rscript assign_asv_taxonomy.R <seqtab_rds> <output_dir> <amplicon> <taxonomy_db> <threads> <project_name> <db_name>

library(dada2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
  stop("Usage: Rscript assign_asv_taxonomy.R <seqtab_rds> <output_dir> <amplicon> <taxonomy_db> <threads> <project_name> <db_name>")
}

seqtab_path <- args[1]
output_dir <- args[2]
amplicon <- args[3]
taxonomy_db <- args[4]
threads <- as.numeric(args[5])
project_name <- args[6]
db_name <- args[7]

cat("\nTAXONOMY ASSIGNMENT\n")
cat("===================\n")
cat("Database:", taxonomy_db, "\n")
cat("Sequence table:", seqtab_path, "\n")

seqtab.nochim <- readRDS(seqtab_path)
cat("ASVs:", ncol(seqtab.nochim), "\n")

# Determine taxonomy levels from the database choice alone (not the amplicon):
# PR2 uses an extended 9-level hierarchy, all others use standard 7-level Kingdom → Species
# (keep in sync with assign_otu_taxonomy.R)
use_extended_levels <- db_name == "PR2"
tax_levels_standard <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_levels_extended <- c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")
tax_levels <- if (use_extended_levels) tax_levels_extended else tax_levels_standard

# Taxonomy parameters: amplicon-specific flags (tryRC, minBoot)
use_tryRC <- amplicon %in% c("18S-AMF", "18S-V4")
use_minBoot <- 50

taxa <- assignTaxonomy(seqtab.nochim, taxonomy_db,
  multithread = threads,
  outputBootstraps = TRUE,
  taxLevels = tax_levels,
  tryRC = use_tryRC,
  minBoot = use_minBoot,
  verbose = TRUE
)

cat("Taxonomy assignment complete\n")

# Save taxonomy RDS objects
taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = file.path(output_dir, "taxID.rds"))
saveRDS(bootout, file = file.path(output_dir, "taxID_bootstrap.rds"))

# Rename bootstrap columns to include "_bootstrap" suffix
colnames(bootout) <- paste0(colnames(bootout), "_bootstrap")

# Build output tables with ASV as a named column (not rowname)
asv_seqs <- rownames(taxout)

both_boot <- data.frame(
  ASV = asv_seqs,
  as.data.frame(t(seqtab.nochim)),
  as.data.frame(taxout),
  as.data.frame(bootout),
  check.names = FALSE
)

both_notax <- data.frame(
  ASV = asv_seqs,
  as.data.frame(t(seqtab.nochim)),
  as.data.frame(taxout),
  check.names = FALSE
)

# Write output files with new naming convention
out_file_taxa <- file.path(output_dir, paste0(project_name, "__combined_sequences_ASVtaxa_", db_name, ".txt"))
out_file_boot <- file.path(output_dir, paste0(project_name, "__combined_sequences_ASVtaxa_bootstrap_", db_name, ".txt"))

write.table(both_notax,
  file = out_file_taxa,
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(both_boot,
  file = out_file_boot,
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nDone.\n")
