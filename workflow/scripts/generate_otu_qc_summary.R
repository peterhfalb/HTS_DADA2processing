#!/usr/bin/env Rscript
# generate_otu_qc_summary.R
# Generates comprehensive OTU QC summary by reading intermediate files from each step
#
# Usage: Rscript generate_otu_qc_summary.R <output_dir> <project_name> <run_itsx>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_otu_qc_summary.R <output_dir> <project_name> <run_itsx>")
}

OUTPUT_DIR  <- args[1]
PROJECT     <- args[2]
# Parse boolean: handle "1", "true", "True", "yes", etc.
RUN_ITSX    <- tolower(args[3]) %in% c("1", "true", "yes")

cat("Generating OTU QC Summary...\n")
cat("Output dir: ", OUTPUT_DIR, "\n")
cat("Project:    ", PROJECT, "\n")
cat("Run ITSx:   ", RUN_ITSX, "\n\n")

# ============================================================================
# Collect counts from each step
# ============================================================================

otu_dir <- file.path(OUTPUT_DIR, "05_asv2otu")

# 1. Input ASVs (from prepare_otu_input)
centroid_input <- file.path(otu_dir, "01_input/Centroid.fasta")
if (file.exists(centroid_input)) {
  input_asvs <- system(paste("grep -c '^>'", centroid_input), intern = TRUE)
  input_asvs <- as.numeric(input_asvs)
  cat("Input ASVs: ", input_asvs, "\n")
} else {
  cat("ERROR: Centroid.fasta not found\n")
  quit(status = 1)
}

# 2. After ITSx (if applicable)
post_itsx <- NA
if (RUN_ITSX) {
  itsx_fasta <- file.path(otu_dir, "02_itsx", paste0("Centroid.ITSx.*.filtered.fasta"))
  # Find the actual file (might be ITS1 or ITS2)
  itsx_files <- Sys.glob(itsx_fasta)
  if (length(itsx_files) > 0) {
    post_itsx <- system(paste("grep -c '^>'", itsx_files[1]), intern = TRUE)
    post_itsx <- as.numeric(post_itsx)
    cat("After ITSx: ", post_itsx, "\n")
  }
}

# 3-5. Read counts from vsearch_cluster rule (clustering, chimera removal, percentages)
vsearch_dir <- file.path(otu_dir, "03_vsearch")
vsearch_counts_file <- file.path(vsearch_dir, paste0(PROJECT, "_vsearch_counts.txt"))

post_cluster <- NA
chimeras_removed <- NA
post_chimera <- NA
pct_chimeric <- NA
pct_asvs_removed <- NA

if (file.exists(vsearch_counts_file)) {
  # Read key=value format file
  counts <- readLines(vsearch_counts_file)
  for (line in counts) {
    parts <- strsplit(line, "=")[[1]]
    if (length(parts) == 2) {
      key <- trimws(parts[1])
      value <- trimws(parts[2])
      if (key == "POST_CLUSTER") post_cluster <- as.numeric(value)
      if (key == "POST_CHIMERA") post_chimera <- as.numeric(value)
      if (key == "CHIMERAS_REMOVED") chimeras_removed <- as.numeric(value)
      if (key == "PCT_CHIMERIC") pct_chimeric <- as.numeric(value)
      if (key == "PCT_ASVS_REMOVED") pct_asvs_removed <- as.numeric(value)
    }
  }
  if (!is.na(post_cluster)) cat("After clustering: ", post_cluster, "\n")
  if (!is.na(chimeras_removed)) cat("Chimeras removed: ", chimeras_removed, "\n")
  if (!is.na(post_chimera)) cat("After chimera removal: ", post_chimera, "\n")
  if (!is.na(pct_chimeric)) cat("% Chimeric: ", pct_chimeric, "%\n")
  if (!is.na(pct_asvs_removed)) cat("% ASVs removed by clustering: ", pct_asvs_removed, "%\n")
} else {
  cat("Warning: VSEARCH counts file not found\n")
}

# 5. Before mumu curation
otutable_path <- file.path(vsearch_dir, paste0(PROJECT, ".otutable"))
if (file.exists(otutable_path)) {
  otus_pre_mumu <- system(paste("tail -n +2", otutable_path, "| wc -l"), intern = TRUE)
  otus_pre_mumu <- as.numeric(otus_pre_mumu)
  cat("OTUs before mumu: ", otus_pre_mumu, "\n")
} else {
  cat("Warning: OTU table not found\n")
  otus_pre_mumu <- NA
}

# 6. After mumu curation
mumu_dir <- file.path(otu_dir, "04_mumu")
mumu_table <- file.path(mumu_dir, paste0(PROJECT, "_mumu_curated.txt"))
otus_post_mumu <- NA
otus_removed_mumu <- NA
pct_removed_mumu <- NA

if (file.exists(mumu_table)) {
  otus_post_mumu <- system(paste("tail -n +2", mumu_table, "| wc -l"), intern = TRUE)
  otus_post_mumu <- as.numeric(otus_post_mumu)
  cat("OTUs after mumu: ", otus_post_mumu, "\n")
} else {
  cat("Warning: mumu curated table not found\n")
}

# OTUs removed by mumu and percentage
if (!is.na(otus_pre_mumu) && !is.na(otus_post_mumu)) {
  otus_removed_mumu <- otus_pre_mumu - otus_post_mumu
  pct_removed_mumu <- round((otus_removed_mumu / otus_pre_mumu) * 100, 1)
  cat("OTUs removed by mumu: ", otus_removed_mumu, " (", pct_removed_mumu, "%)\n")
}

# 7. After taxonomy assignment
tax_dir <- file.path(otu_dir, "05_taxonomy")
tax_files <- list.files(tax_dir, pattern = "_combined.txt$")
tax_file <- if (length(tax_files) > 0) file.path(tax_dir, tax_files[1]) else NA
otus_with_tax <- NA
pct_with_tax <- NA

if (!is.na(tax_file) && file.exists(tax_file)) {
  otus_with_tax <- system(paste("tail -n +2", tax_file, "| wc -l"), intern = TRUE)
  otus_with_tax <- as.numeric(otus_with_tax)
  if (!is.na(otus_post_mumu)) {
    pct_with_tax <- round((otus_with_tax / otus_post_mumu) * 100, 1)
    cat("OTUs with taxonomy: ", otus_with_tax, " (", pct_with_tax, "%)\n")
  } else {
    cat("OTUs with taxonomy: ", otus_with_tax, "\n")
  }
} else {
  cat("Warning: taxonomy file not found\n")
}

# ============================================================================
# Write OTU QC summary
# ============================================================================

output_file <- file.path(otu_dir, "otu_qc_summary.txt")

summary_lines <- c(
  "# ASV to OTU QC Summary",
  paste("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Metric\tValue",
  paste("Input ASVs", input_asvs, sep = "\t")
)

if (!is.na(post_itsx)) {
  summary_lines <- c(summary_lines, paste("After ITSx extraction", post_itsx, sep = "\t"))
}

if (!is.na(post_cluster)) {
  summary_lines <- c(summary_lines, paste("OTUs after clustering (0.97 identity)", post_cluster, sep = "\t"))
}

if (!is.na(pct_asvs_removed)) {
  summary_lines <- c(summary_lines, paste("% ASVs removed by clustering", paste0(pct_asvs_removed, "%"), sep = "\t"))
}

if (!is.na(chimeras_removed)) {
  summary_lines <- c(summary_lines, paste("Chimeras removed", chimeras_removed, sep = "\t"))
}

if (!is.na(pct_chimeric)) {
  summary_lines <- c(summary_lines, paste("% Chimeric sequences", paste0(pct_chimeric, "%"), sep = "\t"))
}

if (!is.na(post_chimera)) {
  summary_lines <- c(summary_lines, paste("OTUs after chimera removal", post_chimera, sep = "\t"))
}

if (!is.na(otus_removed_mumu)) {
  summary_lines <- c(summary_lines, paste("OTUs removed by mumu", otus_removed_mumu, sep = "\t"))
}

if (!is.na(pct_removed_mumu)) {
  summary_lines <- c(summary_lines, paste("% OTUs removed by mumu", paste0(pct_removed_mumu, "%"), sep = "\t"))
}

if (!is.na(otus_post_mumu)) {
  summary_lines <- c(summary_lines, paste("Final OTU count (post-mumu)", otus_post_mumu, sep = "\t"))
}

if (!is.na(otus_with_tax)) {
  if (!is.na(pct_with_tax)) {
    summary_lines <- c(summary_lines, paste("OTUs with taxonomy assigned", paste0(otus_with_tax, " (", pct_with_tax, "%)"), sep = "\t"))
  } else {
    summary_lines <- c(summary_lines, paste("OTUs with taxonomy assigned", otus_with_tax, sep = "\t"))
  }
}

writeLines(summary_lines, output_file)
cat("\n✓ OTU QC summary written to: ", output_file, "\n")
