#!/usr/bin/env Rscript
# DADA2 denoising and taxonomy assignment
# Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads> <project_name> <db_name> <platform> <fwd_reads_only>

library(dada2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads> <project_name> <db_name> <platform> <fwd_reads_only>")
}

output_dir <- args[1]
amplicon <- args[2]
quality <- args[3]
taxonomy_db <- args[4]
threads <- as.numeric(args[5])
project_name <- args[6]
db_name <- args[7]
platform <- tolower(args[8])
# Parse boolean: handle "1", "true", "True", "yes", etc.
fwd_reads_only <- tolower(args[9]) %in% c("1", "true", "yes")

cat("DADA2 Pipeline\n")
cat("==============\n")
cat("Output dir:    ", output_dir, "\n")
cat("Amplicon type: ", amplicon, "\n")
cat("Quality:       ", quality, "\n")
cat("Platform:      ", platform, "\n")
cat("Taxonomy DB:   ", taxonomy_db, "\n")
cat("Project name:  ", project_name, "\n")
cat("DB name:       ", db_name, "\n")
cat("Threads:       ", threads, "\n\n")

# Construct path to filtered reads location (02_primer_trimmed/)
# output_dir is /path/to/FAB2_16S/03_dada2, so parent is /path/to/FAB2_16S
parent_dir <- dirname(output_dir)
filtered_dir <- file.path(parent_dir, "02_primer_trimmed")

cat("Filtered reads directory: ", filtered_dir, "\n\n")

# List files using absolute paths
cat("Finding fastq files...\n")
fnFs <- sort(list.files(filtered_dir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(filtered_dir, pattern = "_R2_001.fastq.gz", full.names = TRUE))
cat("Found ", length(fnFs), " forward and ", length(fnRs), " reverse reads\n")

if (length(fnFs) == 0) {
  stop("No R1 fastq files found")
}
if (!fwd_reads_only && length(fnRs) == 0) {
  stop("No R2 fastq files found")
}

# Extract sample names (strip R1_001 suffix but preserve sample name underscores)
sample.names <- sub("_R1_001\\.fastq\\.gz$", "", basename(fnFs))
cat("Samples: ", paste(sample.names, collapse = ", "), "\n\n")

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
filtered_subdir <- file.path(output_dir, "filtered")
dir.create(filtered_subdir, showWarnings = FALSE, recursive = TRUE)

# Place filtered files
filtFs <- file.path(filtered_subdir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filtered_subdir, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# ==============================================================================
# FILTERING - parameters depend on amplicon type and quality
# ==============================================================================
cat("FILTERING\n")
cat("=========\n")

# Helper function to get filterAndTrim parameters based on amplicon, quality, and platform
get_filter_params <- function(amplicon, quality, platform) {
  # Base parameters (always included)
  base_params <- list(
    maxN = 0,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = min(threads, 4),
    verbose = TRUE
  )

  # Amplicon-specific parameters
  if (amplicon == "16S-V4") {
    if (platform == "aviti") {
      params <- c(base_params, list(truncLen = c(240, 220), truncQ = 5))
      if (quality == "bad") params$maxEE <- c(4, 6) else params$maxEE <- c(2, 2)
      params$minLen <- 100
    } else {  # Illumina
      params <- c(base_params, list(truncLen = c(240, 200), truncQ = 2))
      if (quality == "bad") params$maxEE <- c(4, 6) else params$maxEE <- c(2, 4)
      params$minLen <- 100
    }
  } else if (amplicon == "ITS1" || amplicon == "ITS2") {
    params <- c(base_params, list(truncQ = 2, minLen = 50))
    if (quality == "bad") {
      params$maxEE <- c(4, 4)
      params$truncLen <- c(240, 175)
    } else {
      params$maxEE <- c(2, 2)
    }
  } else if (amplicon == "18S-AMF" || amplicon == "18S-V4") {
    params <- c(base_params, list(truncQ = 5, minLen = 100))
    if (quality == "bad") params$maxEE <- c(4, 6) else params$maxEE <- c(2, 4)
    params$matchIDs <- TRUE
  } else {
    stop("Unknown amplicon type: ", amplicon)
  }

  return(params)
}

# Get parameters and execute filterAndTrim with minimal code duplication
filter_params <- get_filter_params(amplicon, quality, platform)

if (fwd_reads_only) {
  # Forward-only: filter R1 only
  # Adjust paired-end parameters to single-end values (use first element)
  if ("truncLen" %in% names(filter_params) && length(filter_params$truncLen) > 1) {
    filter_params$truncLen <- filter_params$truncLen[1]
  }
  if ("maxEE" %in% names(filter_params) && length(filter_params$maxEE) > 1) {
    filter_params$maxEE <- filter_params$maxEE[1]
  }
  out <- do.call(filterAndTrim, c(list(fnFs, filtFs), filter_params))
} else {
  # Paired-end: filter R1 and R2
  out <- do.call(filterAndTrim, c(list(fnFs, filtFs, fnRs, filtRs), filter_params))
}

head(out)
cat("\n")

# ==============================================================================
# DEREPLICATION
# ==============================================================================
cat("DEREPLICATION\n")
cat("=============\n")

derep_forward <- derepFastq(filtFs, verbose = TRUE)
names(derep_forward) <- sample.names

# ==============================================================================
# ERROR LEARNING & DENOISING
# ==============================================================================
cat("ERROR LEARNING\n")
cat("==============\n")

errF <- learnErrors(derep_forward, multithread = threads, randomize = TRUE)

# Save error models for QC figures
saveRDS(errF, file.path(output_dir, "errF.rds"))

cat("DADA2 INFERENCE\n")
cat("===============\n")

dadaFs <- dada(derep_forward, err = errF, multithread = threads, pool = "pseudo")

# Process reverse reads only if not in forward-only mode
if (!fwd_reads_only) {
  derep_reverse <- derepFastq(filtRs, verbose = TRUE)
  names(derep_reverse) <- sample.names

  errR <- learnErrors(derep_reverse, multithread = threads, randomize = TRUE)
  saveRDS(errR, file.path(output_dir, "errR.rds"))

  dadaRs <- dada(derep_reverse, err = errR, multithread = threads, pool = "pseudo")
}

# ==============================================================================
# MERGING - parameters depend on amplicon type
# ==============================================================================
if (!fwd_reads_only) {
  cat("MERGING PAIRED READS\n")
  cat("====================\n")

  if (amplicon == "16S-V4") {
    merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
      trimOverhang = TRUE, minOverlap = 20
    )
  } else if (amplicon == "ITS1" || amplicon == "ITS2") {
    merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
      trimOverhang = TRUE, minOverlap = 10
    )
  } else if (amplicon == "18S-AMF") {
    merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
      trimOverhang = TRUE, minOverlap = 10
    )
  } else if (amplicon == "18S-V4") {
    merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
      trimOverhang = TRUE, minOverlap = 10
    )
  }
} else {
  # Forward-only mode: skip merging, use forward reads directly
  cat("FORWARD-ONLY MODE\n")
  cat("=================\n")
  errR <- NULL
  saveRDS(errR, file.path(output_dir, "errR.rds"))
  cat("Reverse reads skipped. Sequence table will be built from forward reads only.\n")
}

# ==============================================================================
# SEQUENCE TABLE & CHIMERA REMOVAL
# ==============================================================================
cat("SEQUENCE TABLE CONSTRUCTION\n")
cat("===========================\n")

if (fwd_reads_only) {
  seqtab <- makeSequenceTable(dadaFs)
  cat("Built from forward reads only.\n")
} else {
  seqtab <- makeSequenceTable(merged_amplicons)
}
cat("Dimensions:", dim(seqtab), "\n")
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))

cat("\nCHIMERA REMOVAL\n")
cat("===============\n")

if (platform == "aviti") {
  # Aviti-specific chimera removal: per-sample method with higher abundance threshold
  # This addresses Aviti's overabundance of high-quality error-induced reads
  # See: Gould et al. preprint on Aviti processing best practices
  cat("Using Aviti chimera detection settings:\n")
  cat("  Method: per-sample\n")
  cat("  minFoldParentOverAbundance: 8 (vs default 2)\n\n")
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = "per-sample",
    minFoldParentOverAbundance = 8,
    multithread = threads,
    verbose = TRUE
  )
} else {
  # Illumina (default) chimera removal: consensus method with default abundance threshold
  cat("Using Illumina chimera detection settings:\n")
  cat("  Method: consensus\n")
  cat("  minFoldParentOverAbundance: 2 (default)\n\n")
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = "consensus",
    multithread = threads,
    verbose = TRUE
  )
}

cat("Dimensions after chimera removal:", dim(seqtab.nochim), "\n")
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

# Export sequences
cat("\nEXPORTING SEQUENCES\n")
cat("===================\n")
uniquesToFasta(seqtab.nochim, fout = file.path(output_dir, "sequences.fasta"))
cat("Sequences exported to sequences.fasta\n")

# ==============================================================================
# SUMMARY STATISTICS
# ==============================================================================
cat("\nSUMMARY STATISTICS\n")
cat("==================\n")

getN <- function(x) sum(getUniques(x))

summary_tab <- data.frame(
  row.names = sample.names,
  dada2_input = out[, 1],
  filtered = out[, 2],
  dada_f = sapply(dadaFs, getN),
  dada_r = if (fwd_reads_only) NA else sapply(dadaRs, getN),
  merged = if (fwd_reads_only) NA else sapply(merged_amplicons, getN),
  nonchim = rowSums(seqtab.nochim)
)

print(summary_tab)
write.table(summary_tab,
  file = file.path(output_dir, "sequence_process_summary.txt"),
  sep = "\t", quote = FALSE
)

# ==============================================================================
# TAXONOMY ASSIGNMENT
# ==============================================================================
cat("\nTAXONOMY ASSIGNMENT\n")
cat("===================\n")
cat("Database:", taxonomy_db, "\n")

seqtab.nochim <- readRDS(file.path(output_dir, "seqtab_nochim.rds"))

# Determine taxonomy levels: only PR2 uses extended hierarchy
# All others use standard 7-level Kingdom → Species
use_extended_levels <- grepl("pr2", taxonomy_db, ignore.case = TRUE)
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
