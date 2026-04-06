#!/usr/bin/env Rscript
# DADA2 denoising and taxonomy assignment
# Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads>

library(dada2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_dada2.R <output_dir> <amplicon> <quality> <taxonomy_db> <threads>")
}

output_dir <- args[1]
amplicon <- args[2]
quality <- args[3]
taxonomy_db <- args[4]
threads <- as.numeric(args[5])

cat("DADA2 Pipeline\n")
cat("==============\n")
cat("Output dir:    ", output_dir, "\n")
cat("Amplicon type: ", amplicon, "\n")
cat("Quality:       ", quality, "\n")
cat("Taxonomy DB:   ", taxonomy_db, "\n")
cat("Threads:       ", threads, "\n\n")

# Construct path to filtered reads location (02_primer_trimmed/)
filtered_dir <- dirname(output_dir)
filtered_dir <- file.path(dirname(filtered_dir), "02_primer_trimmed")

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
    trimOverhang = TRUE, minOverlap = 50
  )
} else if (amplicon == "18S-AMF") {
  merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
    trimOverhang = TRUE, minOverlap = 10, justConcatenate = TRUE
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

# Length filtering for 18S amplicons only
if (amplicon == "18S-AMF" || amplicon == "18S-V4") {
  cat("\nLENGTH FILTERING (18S only)\n")
  cat("===========================\n")
  MINLEN <- 400
  MAXLEN <- 600
  seqlens <- nchar(colnames(seqtab))
  cat("Filtering sequences between", MINLEN, "and", MAXLEN, "bp\n")
  seqtab.filt <- seqtab[, seqlens >= MINLEN & seqlens <= MAXLEN]
  cat("Dimensions after length filter:", dim(seqtab.filt), "\n")
  seqtab <- seqtab.filt
}

cat("\nCHIMERA REMOVAL\n")
cat("===============\n")

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
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

# Taxonomy parameters depend on amplicon type
if (amplicon == "16S-V4") {
  taxa <- assignTaxonomy(seqtab.nochim, taxonomy_db,
    multithread = TRUE, outputBootstraps = TRUE
  )
} else if (amplicon == "ITS1" || amplicon == "ITS2") {
  taxa <- assignTaxonomy(seqtab.nochim, taxonomy_db,
    multithread = TRUE, outputBootstraps = TRUE
  )
} else if (amplicon == "18S-AMF") {
  taxa <- assignTaxonomy(seqtab.nochim, taxonomy_db,
    tryRC = TRUE,
    multithread = TRUE, outputBootstraps = TRUE
  )
} else if (amplicon == "18S-V4") {
  taxa <- assignTaxonomy(seqtab.nochim, taxonomy_db,
    multithread = TRUE,
    minBoot = 95,
    verbose = TRUE,
    taxLevels = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"),
    outputBootstraps = TRUE
  )
}

cat("Taxonomy assignment complete\n")

# Save taxonomy
taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = file.path(output_dir, "taxID.rds"))
saveRDS(bootout, file = file.path(output_dir, "taxID_bootstrap.rds"))

# Combine sequences with taxonomy
both_boot <- cbind(t(seqtab.nochim), taxa$tax, taxa$boot)
both_notax <- cbind(t(seqtab.nochim), taxa$tax)

write.table(both_boot,
  file = file.path(output_dir, paste0(amplicon, "_combined_sequences_taxa_bootstrap.txt")),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(both_notax,
  file = file.path(output_dir, paste0(amplicon, "_combined_sequences_taxa.txt")),
  sep = "\t", quote = FALSE, col.names = NA
)

cat("\nDone.\n")
