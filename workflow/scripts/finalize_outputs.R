#!/usr/bin/env Rscript
# finalize_outputs.R
# Finalizes the pipeline outputs:
# 1. Copies ASV table to main directory
# 2. Copies OTU table to main directory (if it exists)
# 3. Combines QC tables and writes to main directory
# 4. Generates README for main directory
#
# Usage: Rscript finalize_outputs.R <output_dir> <project_name> <db_name> <skip_otu>

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript finalize_outputs.R <output_dir> <project_name> <db_name> <skip_otu>")
}

OUTPUT_DIR  <- args[1]
PROJECT     <- args[2]
DB_NAME     <- args[3]
# Handle both Snakemake Python booleans ("True"/"False") and shell numbers ("0"/"1")
skip_otu_str <- tolower(args[4])
SKIP_OTU    <- skip_otu_str %in% c("true", "1", "yes")

setwd(OUTPUT_DIR)

cat("Finalizing outputs...\n")
cat("Output dir:  ", OUTPUT_DIR, "\n")
cat("Project:     ", PROJECT, "\n")
cat("DB:          ", DB_NAME, "\n")
cat("Skip OTU:    ", SKIP_OTU, "\n\n")

# ============================================================================
# 1. Copy ASV table with bootstrap values to main directory
# ============================================================================
asv_source <- paste0("03_dada2/", PROJECT, "__combined_sequences_ASVtaxa_bootstrap_", DB_NAME, ".txt")
asv_dest <- paste0(PROJECT, "__combined_sequences_ASVtaxa_bootstrap_", DB_NAME, ".txt")

if (file.exists(asv_source)) {
  file.copy(asv_source, asv_dest, overwrite = TRUE)
  cat("✓ Copied ASV table with bootstrap values to main directory: ", asv_dest, "\n")
} else {
  warning("ASV table not found: ", asv_source)
}

# ============================================================================
# 2. Copy OTU table to main directory (if it exists)
# ============================================================================
otu_source <- paste0("05_asv2otu/", PROJECT, "__OTUs_with_taxonomy_", DB_NAME, ".txt")
otu_dest <- paste0(PROJECT, "__OTUs_with_taxonomy_", DB_NAME, ".txt")

otu_exists <- FALSE
if (!SKIP_OTU && file.exists(otu_source)) {
  file.copy(otu_source, otu_dest, overwrite = TRUE)
  cat("✓ Copied OTU table to main directory: ", otu_dest, "\n")
  otu_exists <- TRUE
} else if (SKIP_OTU) {
  cat("✓ OTU pipeline skipped; not copying OTU table\n")
} else if (!file.exists(otu_source)) {
  warning("OTU table not found (expected if OTU pipeline skipped): ", otu_source)
}

# ============================================================================
# 3. Combine QC tables and write to main directory
# ============================================================================
cat("\nProcessing QC table...\n")

# Read per-sample QC table from DADA2/cutadapt step
qc_per_sample_path <- "04_dada2_QCsummary/qc_summary.txt"
if (file.exists(qc_per_sample_path)) {
  qc_per_sample <- read_tsv(qc_per_sample_path, comment = "#", show_col_types = FALSE)
  cat("✓ Loaded per-sample QC table (", nrow(qc_per_sample), " rows)\n")
} else {
  warning("Per-sample QC table not found: ", qc_per_sample_path)
  qc_per_sample <- data.frame()
}

# Prepare output QC table (start with per-sample)
qc_final <- qc_per_sample

# Append OTU QC summary if OTU was run
otu_qc_path <- "05_asv2otu/otu_qc_summary.txt"
if (!SKIP_OTU && file.exists(otu_qc_path)) {
  cat("✓ Found OTU QC summary\n")

  # Append blank line, section header, and OTU QC section
  blank_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(qc_final)))
  names(blank_row) <- names(qc_final)
  qc_final <- bind_rows(qc_final, blank_row)

  # Add section header
  header_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(qc_final)))
  names(header_row) <- names(qc_final)
  header_row[[1]] <- "# ASV to OTU QC Summary"
  qc_final <- bind_rows(qc_final, header_row)

  # Read OTU QC summary (Metric\tValue format, skip 2 comments + 1 blank line)
  otu_qc <- read_tsv(otu_qc_path, skip = 3, show_col_types = FALSE)

  # Add each OTU metric as a data row
  for (i in seq_len(nrow(otu_qc))) {
    metric_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(qc_final)))
    names(metric_row) <- names(qc_final)
    metric_row[[1]] <- paste(otu_qc[[i, 1]], otu_qc[[i, 2]], sep = ": ")
    qc_final <- bind_rows(qc_final, metric_row)
  }
} else if (!SKIP_OTU) {
  cat("⚠  OTU QC summary file not found (non-critical if OTU pipeline is still running)\n")
}

# Write combined QC table
write.table(qc_final,
  file = "qc_summary.txt",
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
cat("✓ Wrote combined QC table: qc_summary.txt\n")

# ============================================================================
# 4. Generate README for main directory
# ============================================================================
cat("\nGenerating main directory README...\n")

readme_content <- paste0("# ", PROJECT, " — Final Analysis Results\n\n")

readme_content <- paste0(readme_content, "## Summary Output Files\n\n")

readme_content <- paste0(readme_content,
  "### ASV Table (Amplicon Sequence Variants)\n",
  "**File:** `", asv_dest, "`\n\n",
  "Denoised and dereplicated sequences from the DADA2 pipeline.\n\n",
  "**Columns:**\n",
  "- `ASV`: The DNA sequence itself\n",
  "- Sample columns: Read counts per ASV per sample\n",
  "- Taxonomy columns: Kingdom, Phylum, Class, Order, Family, Genus, Species\n",
  "  (or Domain, Supergroup, ... for PR2 database)\n",
  "- Bootstrap columns: Confidence values (0-100) for each taxonomy rank\n\n",
  "**Interpretation:**\n",
  "- Each row is a unique sequence variant\n",
  "- Use for: diversity analysis, relative abundance plots, comparative studies\n\n"
)

if (otu_exists) {
  readme_content <- paste0(readme_content,
    "### OTU Table (Operational Taxonomic Units)\n",
    "**File:** `", otu_dest, "`\n\n",
    "Clustered and curated OTUs from VSEARCH (97% identity) and mumu curation.\n\n",
    "**Columns:**\n",
    "- `OTU_ID`: OTU identifier\n",
    "- Sample columns: Read counts per OTU per sample\n",
    "- Taxonomy columns: Kingdom through Species (matching ASV table taxonomy)\n",
    "- Bootstrap columns: Confidence values for taxonomy assignments\n\n",
    "**Interpretation:**\n",
    "- Each row is an OTU (97% similar sequences grouped together)\n",
    "- Represents operational taxa rather than true biological species\n",
    "- Use for: traditional OTU-based analyses, clustering, rarefaction curves\n",
    "- **Note:** OTUs may contain sequences from similar but distinct taxa\n\n"
  )
}

readme_content <- paste0(readme_content,
  "### QC Summary Table\n",
  "**File:** `qc_summary.txt`\n\n",
  "Quality control metrics from all pipeline steps.\n\n"
)

if (!SKIP_OTU) {
  readme_content <- paste0(readme_content,
    "**Section 1: Per-Sample Metrics**\n",
    "For each sample:\n",
    "- Read counts entering each step (cutadapt, DADA2, merging, chimera removal)\n",
    "- Percentage of reads passing through each filter\n",
    "- Quality metrics from DADA2 inference\n\n",
    "**Section 2: OTU Processing Summary**\n",
    "Pipeline-wide OTU clustering statistics:\n",
    "- `input_asvs`: Unique ASVs entering OTU step\n",
    "- `post_itsx`: ASVs after ITSx extraction (ITS datasets only)\n",
    "- `post_cluster`: Sequences after VSEARCH 97% clustering\n",
    "- `post_chimera`: Sequences after chimera removal\n",
    "- `post_mumu`: OTUs after mumu curation\n",
    "- `classified_otus`: OTUs with successful taxonomy assignment\n\n"
  )
} else {
  readme_content <- paste0(readme_content,
    "Covers per-sample metrics from cutadapt and DADA2 processing steps.\n\n"
  )
}

readme_content <- paste0(readme_content,
  "## Additional Output Directories\n\n",
  "- `03_dada2/`: Full DADA2 denoising pipeline outputs (error models, sequence tables, taxonomy RDS files)\n",
  "- `04_dada2_QCsummary/`: Quality control figures and statistics (quality profiles, error models, read tracking)\n"
)

if (!SKIP_OTU) {
  readme_content <- paste0(readme_content,
    "- `05_asv2otu/`: OTU clustering workflow details (VSEARCH centroids, mumu curation logs, intermediate tables)\n"
  )
}

readme_content <- paste0(readme_content,
  "\n## Choosing Between ASV and OTU Tables\n\n",
  "**Use ASV table if:**\n",
  "- You want to preserve exact sequence variants\n",
  "- You're interested in rare sequence variants\n",
  "- You're performing strain-level analysis\n\n",
  "**Use OTU table if:**\n",
  "- You need traditional OTU-based analyses\n",
  "- You want to reduce computational complexity\n",
  "- You're comparing with older studies using OTU definitions\n\n"
)

writeLines(readme_content, "README.txt")
cat("✓ Generated README.txt\n")

cat("\n✓ Pipeline finalization complete!\n")
