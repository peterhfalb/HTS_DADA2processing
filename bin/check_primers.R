#!/usr/bin/env Rscript
# Primer validation: check if primers found in raw FASTQ match expected primers
# Usage: Rscript check_primers.R <fastq_dir> <expected_fwd_primer> <expected_rev_primer> <amplicon>

# Primer reference database
PRIMER_DB <- list(
  # ITS1
  "ITS1F" = list(name = "ITS1F", region = "ITS1", direction = "Forward", seq = "CTTGGTCATTTAGAGGAAGTAA"),
  "ITS2" = list(name = "ITS2", region = "ITS1", direction = "Reverse", seq = "GCTGCGTTCTTCATCGATGC"),
  "ITS1F_KYO1" = list(name = "ITS1F_KYO1", region = "ITS1", direction = "Forward", seq = "CTHGGTCATTTAGAGGAASTAA"),
  "ITS2_KYO2" = list(name = "ITS2_KYO2", region = "ITS1", direction = "Reverse", seq = "TTYRCTRYGTTCTTCATC"),

  # ITS2
  "5.8S-FUN" = list(name = "5.8S-FUN", region = "ITS2", direction = "Forward", seq = "AACTTTYRRCAAYGGATCWCT"),
  "ITS4-Fun" = list(name = "ITS4-Fun", region = "ITS2", direction = "Reverse", seq = "AGCCTCCGCTTATTGATATGCTTAART"),
  "fITS7_1" = list(name = "fITS7_1", region = "ITS2", direction = "Forward", seq = "GTGARTCATCGAATCTTG"),
  "ITS4_1" = list(name = "ITS4_1", region = "ITS2", direction = "Reverse", seq = "TCCTCCGCTTATTGATATGC"),
  "5.8SR" = list(name = "5.8SR", region = "ITS2", direction = "Forward", seq = "TCGATGAAGAACGCAGCG"),
  "ITS4" = list(name = "ITS4", region = "ITS2", direction = "Reverse", seq = "TCCTCCGCTTATTGATATGC"),

  # 18S SSU (AMF)
  "WANDA_1" = list(name = "WANDA_1", region = "SSU", direction = "Forward", seq = "CAGCCGCGGTAATTCCAGCT"),
  "AML2_1" = list(name = "AML2_1", region = "SSU", direction = "Reverse", seq = "GAACCCAAACACTTTGGTTTCC"),

  # 16S V4
  "515f" = list(name = "515f", region = "16S_V4", direction = "Forward", seq = "GTGYCAGCMGCCGCGGTAA"),
  "806r" = list(name = "806r", region = "16S_V4", direction = "Reverse", seq = "GGACTACNVGGGTWTCTAAT"),
  "27F" = list(name = "27F", region = "16S_V1-V9", direction = "Forward", seq = "AGRTTTGATYMTGGCTCAG"),
  "1492R" = list(name = "1492R", region = "16S_V1-V9", direction = "Reverse", seq = "RGYTACCTTGTTACGACTT"),

  # 18S (Metazoans)
  "F04" = list(name = "F04", region = "18S", direction = "Forward", seq = "GCTTGTCTCAAAGATTAAGCC"),
  "R22" = list(name = "R22", region = "18S", direction = "Reverse", seq = "GCCTGCTGCCTTCCTTGGA"),

  # 18S (Microbial eukaryotes)
  "616f" = list(name = "616f", region = "18S", direction = "Forward", seq = "TTAAARVGYTCGTAGTYG"),
  "1132r" = list(name = "1132r", region = "18S", direction = "Reverse", seq = "CCGTCAATTHCTTTYAART")
)

# IUPAC ambiguity codes mapping
IUPAC_CODES <- list(
  "R" = c("A", "G"),
  "Y" = c("C", "T"),
  "S" = c("G", "C"),
  "W" = c("A", "T"),
  "K" = c("G", "T"),
  "M" = c("A", "C"),
  "B" = c("C", "G", "T"),
  "D" = c("A", "G", "T"),
  "H" = c("A", "C", "T"),
  "V" = c("A", "C", "G"),
  "N" = c("A", "C", "G", "T")
)

# Function to expand IUPAC codes to all possible sequences
expand_iupac <- function(seq) {
  seq_upper <- toupper(seq)
  chars <- strsplit(seq_upper, "")[[1]]

  # Check if any IUPAC codes are present
  if (!any(chars %in% names(IUPAC_CODES))) {
    return(c(seq_upper))
  }

  # Recursively expand IUPAC codes
  result <- c("")
  for (char in chars) {
    if (char %in% names(IUPAC_CODES)) {
      options <- IUPAC_CODES[[char]]
      new_result <- c()
      for (opt in options) {
        new_result <- c(new_result, paste0(result, opt))
      }
      result <- new_result
    } else {
      result <- paste0(result, char)
    }
  }
  return(result)
}

# Function to find primer matches in FASTQ (allowing for some mismatch/IUPAC)
find_primer_in_fastq <- function(fastq_file, primer_seq, max_mismatches = 1, n_reads_to_check = 1000) {
  tryCatch({
    # Read FASTQ file (every 4th line is sequence)
    lines <- readLines(fastq_file)

    # Extract sequence lines (every 4th line starting at line 2)
    n_total <- length(lines) %/% 4
    seq_indices <- seq(2, length(lines), by = 4)
    sequences <- lines[seq_indices[1:min(n_reads_to_check, length(seq_indices))]]

    primer_upper <- toupper(primer_seq)
    primer_len <- nchar(primer_upper)

    # Count exact matches
    exact_count <- sum(grepl(paste0("^", primer_upper), sequences, ignore.case = TRUE))

    # Count matches allowing for IUPAC ambiguity
    # Build regex pattern from IUPAC codes
    pattern <- ""
    for (char in strsplit(primer_upper, "")[[1]]) {
      if (char %in% names(IUPAC_CODES)) {
        options <- IUPAC_CODES[[char]]
        pattern <- paste0(pattern, "[", paste(options, collapse = ""), "]")
      } else {
        pattern <- paste0(pattern, char)
      }
    }

    iupac_count <- sum(grepl(paste0("^", pattern), sequences, ignore.case = TRUE))

    list(
      exact_matches = exact_count,
      iupac_matches = iupac_count,
      total_checked = length(sequences),
      pct_exact = 100 * exact_count / length(sequences),
      pct_iupac = 100 * iupac_count / length(sequences)
    )
  }, error = function(e) {
    list(exact_matches = NA, iupac_matches = NA, total_checked = NA, pct_exact = NA, pct_iupac = NA, error = TRUE)
  })
}

# Main function
validate_primers <- function(fastq_dir, expected_fwd, expected_rev, amplicon) {
  cat("Primer Validation Check\n")
  cat("=======================\n")
  cat("Fastq directory: ", fastq_dir, "\n")
  cat("Amplicon:        ", amplicon, "\n")
  cat("Expected fwd:    ", expected_fwd, "\n")
  cat("Expected rev:    ", expected_rev, "\n\n")

  # Find first R1 and R2 files
  r1_files <- sort(list.files(fastq_dir, pattern = "_R1_001.fastq.gz$", full.names = TRUE))
  r2_files <- sort(list.files(fastq_dir, pattern = "_R2_001.fastq.gz$", full.names = TRUE))

  if (length(r1_files) == 0 || length(r2_files) == 0) {
    cat("ERROR: No R1 or R2 fastq files found in", fastq_dir, "\n")
    return(invisible(NULL))
  }

  # Use first file
  r1_file <- r1_files[1]
  r2_file <- r2_files[1]

  cat("Checking file: ", basename(r1_file), "\n\n")

  # Check forward primer
  cat("FORWARD PRIMER ANALYSIS\n")
  cat("------------------------\n")
  cat("Primer sequence: ", expected_fwd, "\n")

  # Find matching primer in DB
  matched_fwd <- NULL
  for (primer_name in names(PRIMER_DB)) {
    if (toupper(PRIMER_DB[[primer_name]]$seq) == toupper(expected_fwd)) {
      matched_fwd <- PRIMER_DB[[primer_name]]$name
      cat("Matches database primer: ", matched_fwd, "\n")
      break
    }
  }

  if (is.null(matched_fwd)) {
    cat("Note: Primer not found in reference database\n")
  }

  # Search in R1 file
  cat("Searching in R1 reads (first 1000)...\n")
  fwd_results <- find_primer_in_fastq(r1_file, expected_fwd)

  if (!is.null(fwd_results$error)) {
    cat("ERROR: Could not read fastq file\n")
  } else {
    cat("Exact matches:   ", fwd_results$exact_matches, "/", fwd_results$total_checked,
      " (", round(fwd_results$pct_exact, 2), "%)\n")
    cat("IUPAC matches:   ", fwd_results$iupac_matches, "/", fwd_results$total_checked,
      " (", round(fwd_results$pct_iupac, 2), "%)\n")
  }

  # Check reverse primer
  cat("\nREVERSE PRIMER ANALYSIS\n")
  cat("------------------------\n")
  cat("Primer sequence: ", expected_rev, "\n")

  # Find matching primer in DB
  matched_rev <- NULL
  for (primer_name in names(PRIMER_DB)) {
    if (toupper(PRIMER_DB[[primer_name]]$seq) == toupper(expected_rev)) {
      matched_rev <- PRIMER_DB[[primer_name]]$name
      cat("Matches database primer: ", matched_rev, "\n")
      break
    }
  }

  if (is.null(matched_rev)) {
    cat("Note: Primer not found in reference database\n")
  }

  # Search in R2 file
  cat("Searching in R2 reads (first 1000)...\n")
  rev_results <- find_primer_in_fastq(r2_file, expected_rev)

  if (!is.null(rev_results$error)) {
    cat("ERROR: Could not read fastq file\n")
  } else {
    cat("Exact matches:   ", rev_results$exact_matches, "/", rev_results$total_checked,
      " (", round(rev_results$pct_exact, 2), "%)\n")
    cat("IUPAC matches:   ", rev_results$iupac_matches, "/", rev_results$total_checked,
      " (", round(rev_results$pct_iupac, 2), "%)\n")
  }

  # Summary and recommendations
  cat("\nSUMMARY & RECOMMENDATIONS\n")
  cat("===========================\n")

  fwd_ok <- !is.null(fwd_results$error) && fwd_results$pct_iupac > 50
  rev_ok <- !is.null(rev_results$error) && rev_results$pct_iupac > 50

  if (fwd_ok && rev_ok) {
    cat("✓ Both primers detected at high frequency. Looks good!\n")
  } else if (fwd_ok || rev_ok) {
    cat("⚠ One primer detected, other not found or at low frequency.\n")
    cat("  Please verify your primer sequences are correct.\n")
  } else {
    cat("✗ Neither primer detected at significant frequency.\n")
    cat("  CAUTION: Check that:\n")
    cat("    1. Primer sequences match your input files\n")
    cat("    2. Input files are in the correct directory\n")
    cat("    3. Files are in the expected demultiplexed format\n")
  }

  cat("\nDo NOT proceed with SLURM submission until primers are confirmed!\n")
}

# Parse arguments and run
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript check_primers.R <fastq_dir> <expected_fwd_primer> <expected_rev_primer> <amplicon>")
}

fastq_dir <- args[1]
expected_fwd <- args[2]
expected_rev <- args[3]
amplicon <- args[4]

validate_primers(fastq_dir, expected_fwd, expected_rev, amplicon)
