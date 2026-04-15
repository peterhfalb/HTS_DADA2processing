#!/usr/bin/env Rscript
# DADA2 denoising and taxonomy assignment
# Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads> <project_name> <db_name> <platform>

library(dada2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads> <project_name> <db_name> <platform>")
}

output_dir <- args[1]
amplicon <- args[2]
quality <- args[3]
taxonomy_db <- args[4]
threads <- as.numeric(args[5])
project_name <- args[6]
db_name <- args[7]
platform <- tolower(args[8])

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

if (length(fnFs) == 0 || length(fnRs) == 0) {
  stop("No fastq files found")
}

# Extract sample names (ProjectName_SampleName format)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
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

if (amplicon == "16S-V4") {
  if (platform == "aviti") {
    # Aviti platform-specific settings for 16S-V4
    if (quality == "good") {
      out <- filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        truncLen = c(240, 220),
        maxN = 0,
        maxEE = c(2, 2),
        minLen = 100,
        truncQ = 2,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads
      )
    } else if (quality == "bad") {
      out <- filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        truncLen = c(260, 260),
        maxN = 0,
        maxEE = c(4, 6),
        minLen = 100,
        truncQ = 2,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads
      )
    }
  } else {
    # Illumina (default) platform settings for 16S-V4
    if (quality == "good") {
      out <- filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        truncLen = c(240, 200),
        maxN = 0,
        maxEE = c(2, 4),
        minLen = 100,
        truncQ = 2,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads
      )
    } else if (quality == "bad") {
      out <- filterAndTrim(
        fnFs, filtFs, fnRs, filtRs,
        truncLen = c(240, 200),
        maxN = 0,
        maxEE = c(4, 6),
        minLen = 100,
        truncQ = 2,
        rm.phix = TRUE,
        compress = TRUE,
        multithread = threads
      )
    }
  }
} else if (amplicon == "ITS1" || amplicon == "ITS2") {
  if (quality == "good") {
    out <- filterAndTrim(
      fnFs, filtFs, fnRs, filtRs,
      maxN = 0,
      maxEE = c(2, 2),
      minLen = 50,
      truncQ = 2,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = threads
    )
  } else if (quality == "bad") {
    out <- filterAndTrim(
      fnFs, filtFs, fnRs, filtRs,
      maxN = 0,
      maxEE = c(4, 4),
      minLen = 50,
      truncLen=c(240,175),
      truncQ = 2,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = threads
    )
  }
} else if (amplicon == "18S-AMF" || amplicon == "18S-V4") {
  if (quality == "good") {
    out <- filterAndTrim(
      fnFs, filtFs, fnRs, filtRs,
      truncQ = 5,
      minLen = 100,
      maxEE = c(2, 4),
      matchIDs = TRUE,
      maxN = 0,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = threads,
      verbose = TRUE
    )
  } else if (quality == "bad") {
    out <- filterAndTrim(
      fnFs, filtFs, fnRs, filtRs,
      truncQ = 5,
      minLen = 100,
      maxEE = c(4, 6),
      matchIDs = TRUE,
      maxN = 0,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = threads,
      verbose = TRUE
    )
  }
} else {
  stop("Unknown amplicon type: ", amplicon)
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
derep_reverse <- derepFastq(filtRs, verbose = TRUE)
names(derep_reverse) <- sample.names

# ==============================================================================
# ERROR LEARNING & DENOISING
# ==============================================================================
cat("ERROR LEARNING\n")
cat("==============\n")

errF <- learnErrors(derep_forward, multithread = threads, randomize = TRUE)
errR <- learnErrors(derep_reverse, multithread = threads, randomize = TRUE)

# Save error models for QC figures
saveRDS(errF, file.path(output_dir, "errF.rds"))
saveRDS(errR, file.path(output_dir, "errR.rds"))

cat("DADA2 INFERENCE\n")
cat("===============\n")

dadaFs <- dada(derep_forward, err = errF, multithread = threads, pool = "pseudo")
dadaRs <- dada(derep_reverse, err = errR, multithread = threads, pool = "pseudo")

# ==============================================================================
# MERGING - parameters depend on amplicon type
# ==============================================================================
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

# ==============================================================================
# SEQUENCE TABLE & CHIMERA REMOVAL
# ==============================================================================
cat("SEQUENCE TABLE CONSTRUCTION\n")
cat("===========================\n")

seqtab <- makeSequenceTable(merged_amplicons)
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
    multithread = TRUE,
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
    multithread = TRUE,
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
  dada_r = sapply(dadaRs, getN),
  merged = sapply(merged_amplicons, getN),
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
  multithread = TRUE,
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
