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
RUN_ITSX    <- args[3] == "1"

cat("Generating OTU QC Summary...\n")
cat("Output dir: ", OUTPUT_DIR, "\n")
cat("Project:    ", PROJECT, "\n")
cat("Run ITSx:   ", RUN_ITSX, "\n\n")

# ============================================================================
# Collect counts from each step
# ============================================================================

otu_dir <- file.path(OUTPUT_DIR, "05_otu")

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

# 3. After clustering (before chimera check)
vsearch_dir <- file.path(otu_dir, "03_vsearch")
clustered_fasta <- file.path(vsearch_dir, "clustered.fasta")
if (file.exists(clustered_fasta)) {
  # We need to count from the .uc file or the centroids
  # The easiest is to count from the .uc file - each S record is a centroid
  uc_file <- file.path(vsearch_dir, paste0(PROJECT, ".uc"))
  if (file.exists(uc_file)) {
    post_cluster <- system(paste("grep '^S' ", uc_file, "| wc -l"), intern = TRUE)
    post_cluster <- as.numeric(post_cluster)
    cat("After clustering: ", post_cluster, "\n")
  } else {
    cat("Warning: .uc file not found for cluster count\n")
    post_cluster <- NA
  }
} else {
  cat("Warning: clustered.fasta not found\n")
  post_cluster <- NA
}

# 4. After chimera removal
nochim_fasta <- file.path(vsearch_dir, "nochim.fasta")
if (file.exists(nochim_fasta)) {
  post_chimera <- system(paste("grep -c '^>'", nochim_fasta), intern = TRUE)
  post_chimera <- as.numeric(post_chimera)
  cat("After chimera removal: ", post_chimera, "\n")
} else {
  cat("Warning: nochim.fasta not found\n")
  post_chimera <- NA
}

# Chimeras removed
chimeras_removed <- NA
if (!is.na(post_cluster) && !is.na(post_chimera)) {
  chimeras_removed <- post_cluster - post_chimera
  cat("Chimeras removed: ", chimeras_removed, "\n")
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
if (file.exists(mumu_table)) {
  otus_post_mumu <- system(paste("tail -n +2", mumu_table, "| wc -l"), intern = TRUE)
  otus_post_mumu <- as.numeric(otus_post_mumu)
  cat("OTUs after mumu: ", otus_post_mumu, "\n")
} else {
  cat("Warning: mumu curated table not found\n")
  otus_post_mumu <- NA
}

# OTUs removed by mumu
otus_removed_mumu <- NA
if (!is.na(otus_pre_mumu) && !is.na(otus_post_mumu)) {
  otus_removed_mumu <- otus_pre_mumu - otus_post_mumu
  cat("OTUs removed by mumu: ", otus_removed_mumu, "\n")
}

# 7. After taxonomy assignment
tax_dir <- file.path(otu_dir, "05_taxonomy")
tax_file <- file.path(tax_dir, list.files(tax_dir, pattern = "_combined.txt$")[1])
if (!is.na(tax_file) && file.exists(tax_file)) {
  otus_with_tax <- system(paste("tail -n +2", tax_file, "| wc -l"), intern = TRUE)
  otus_with_tax <- as.numeric(otus_with_tax)
  cat("OTUs with taxonomy: ", otus_with_tax, "\n")
} else {
  cat("Warning: taxonomy file not found\n")
  otus_with_tax <- NA
}

# ============================================================================
# Write OTU QC summary
# ============================================================================

output_file <- file.path(otu_dir, "otu_qc_summary.txt")

summary_lines <- c(
  "# OTU QC Summary",
  paste("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Metric\tValue",
  paste("Input ASVs", input_asvs, sep = "\t")
)

if (!is.na(post_itsx)) {
  summary_lines <- c(summary_lines, paste("After ITSx extraction", post_itsx, sep = "\t"))
}

if (!is.na(post_cluster)) {
  summary_lines <- c(summary_lines, paste("After clustering (0.97 identity)", post_cluster, sep = "\t"))
}

if (!is.na(chimeras_removed)) {
  summary_lines <- c(summary_lines, paste("Chimeras removed", chimeras_removed, sep = "\t"))
}

if (!is.na(post_chimera)) {
  summary_lines <- c(summary_lines, paste("OTUs after chimera removal", post_chimera, sep = "\t"))
}

if (!is.na(otus_removed_mumu)) {
  summary_lines <- c(summary_lines, paste("OTUs removed by mumu", otus_removed_mumu, sep = "\t"))
}

if (!is.na(otus_post_mumu)) {
  summary_lines <- c(summary_lines, paste("Final OTU count (post-mumu)", otus_post_mumu, sep = "\t"))
}

if (!is.na(otus_with_tax)) {
  summary_lines <- c(summary_lines, paste("OTUs with taxonomy assigned", otus_with_tax, sep = "\t"))
}

writeLines(summary_lines, output_file)
cat("\n✓ OTU QC summary written to: ", output_file, "\n")
