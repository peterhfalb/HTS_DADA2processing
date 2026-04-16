#!/usr/bin/env Rscript
# QC figures and stats generation for DADA2 pipeline
# Usage: Rscript qc_figures.R <qc_outdir> <amplicon> <platform> <dada2_summary_path> <dada2_dir> <primer_trimmed_dir>

library(dada2)
library(jsonlite)
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript qc_figures.R <qc_outdir> <amplicon> <platform> <dada2_summary_path> <dada2_dir> <primer_trimmed_dir>")
}

qc_outdir <- args[1]
amplicon <- args[2]
platform <- tolower(args[3])
dada2_summary_path <- args[4]
dada2_dir <- args[5]
primer_trimmed_dir <- args[6]

cat("QC Figures Generation\n")
cat("====================\n")
cat("QC output dir:        ", qc_outdir, "\n")
cat("Amplicon type:        ", amplicon, "\n")
cat("Platform:             ", platform, "\n")
cat("DADA2 summary:        ", dada2_summary_path, "\n")
cat("DADA2 dir:            ", dada2_dir, "\n")
cat("Primer-trimmed dir:   ", primer_trimmed_dir, "\n\n")

# Ensure output directories exist
dir.create(file.path(qc_outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PARSE TRIMMOMATIC LOGS (for Aviti adapter trimming)
# ==============================================================================
parse_trimmomatic_log <- function(log_path) {
  tryCatch({
    lines <- readLines(log_path)
    # Trimmomatic log format: "TrimmomaticPE: Processed 1000000 paired reads in 50.23 seconds (19882.15 reads/sec)"
    # Extract pairs processed and various outcomes

    # Look for the main summary line
    summary_line <- grep("TrimmomaticPE:", lines, value = TRUE)

    if (length(summary_line) == 0) {
      return(list(input = NA, output = NA, pct_removed = NA, pct_r1 = NA, pct_r2 = NA))
    }

    summary_line <- summary_line[length(summary_line)]  # Use last match if multiple

    # Extract input reads: "Processed X paired reads"
    input_match <- regexpr("Processed\\s+([0-9]+)\\s+paired", summary_line)
    input_reads <- NA
    if (input_match > 0) {
      input_reads <- as.numeric(sub("^.*Processed\\s+([0-9]+).*", "\\1", summary_line))
    }

    # Parse the per-read statistics to get survival count
    # Lines with format: "   Both surviving: 1000 (100.00%)"
    both_surviving <- grep("Both surviving", lines, value = TRUE)
    output_reads <- NA
    if (length(both_surviving) > 0) {
      both_line <- both_surviving[length(both_surviving)]
      output_match <- regexpr("([0-9]+)\\s+\\(", both_line)
      if (output_match > 0) {
        output_reads <- as.numeric(sub("^.*\\s+([0-9]+)\\s+\\(.*", "\\1", both_line))
      }
    }

    list(
      input = input_reads,
      output = output_reads,
      r1_with_adapter = NA,
      r2_with_adapter = NA,
      pct_r1 = NA,
      pct_r2 = NA,
      pct_removed = ifelse(!is.na(input_reads) && !is.na(output_reads) && input_reads > 0,
                           100 * (input_reads - output_reads) / input_reads, NA)
    )
  }, error = function(e) {
    cat("      ERROR:", e$message, "\n")
    list(input = NA, output = NA, r1_with_adapter = NA, r2_with_adapter = NA, pct_r1 = NA, pct_r2 = NA, pct_removed = NA)
  })
}

# ==============================================================================
# PARSE CUTADAPT JSON FILES
# ==============================================================================
cat("Parsing trimming data...\n")

parse_cutadapt_json <- function(json_path) {
  tryCatch({
    data <- jsonlite::fromJSON(json_path)

    # Extract read counts - handle different cutadapt versions
    if (!is.null(data$read_counts)) {
      rc <- data$read_counts
    } else if (!is.null(data$statistics)) {
      rc <- data$statistics
    } else {
      rc <- list(input = NA, output = NA)
    }

    total_reads <- if (!is.null(rc$input)) rc$input else NA
    output_reads <- if (!is.null(rc$output)) rc$output else NA

    # Try to extract adapter counts from various possible field names
    r1_adapter <- NA
    r2_adapter <- NA

    if (!is.null(rc$read1_with_adapter)) {
      r1_adapter <- rc$read1_with_adapter
    } else if (!is.null(rc$reads_read1_with_trimmed_adapter)) {
      r1_adapter <- rc$reads_read1_with_trimmed_adapter
    }

    if (!is.null(rc$read2_with_adapter)) {
      r2_adapter <- rc$read2_with_adapter
    } else if (!is.null(rc$reads_read2_with_trimmed_adapter)) {
      r2_adapter <- rc$reads_read2_with_trimmed_adapter
    }

    list(
      input = total_reads,
      output = output_reads,
      r1_with_adapter = r1_adapter,
      r2_with_adapter = r2_adapter,
      pct_r1 = ifelse(!is.na(r1_adapter) && !is.na(total_reads) && total_reads > 0, 100 * r1_adapter / total_reads, NA),
      pct_r2 = ifelse(!is.na(r2_adapter) && !is.na(total_reads) && total_reads > 0, 100 * r2_adapter / total_reads, NA),
      pct_removed = ifelse(!is.na(total_reads) && !is.na(output_reads) && total_reads > 0, 100 * (total_reads - output_reads) / total_reads, NA)
    )
  }, error = function(e) {
    cat("      ERROR:", e$message, "\n")
    list(input = NA, output = NA, r1_with_adapter = NA, r2_with_adapter = NA, pct_r1 = NA, pct_r2 = NA, pct_removed = NA)
  })
}

# Parse adapter logs (format depends on platform)
base_output_dir <- dirname(primer_trimmed_dir)
cat("Base output directory:", base_output_dir, "\n")

adapter_log_dir <- if (platform == "aviti") {
  file.path(base_output_dir, ".logs")
} else {
  file.path(base_output_dir, "01_adapter", "01_logs")
}

cat("Looking for adapter logs in:", adapter_log_dir, "\n")
cat("  Directory exists:", dir.exists(adapter_log_dir), "\n")

adapter_stats <- list()

if (platform == "aviti") {
  cat("Platform is Aviti: looking for Trimmomatic stderr logs in .logs/\n")
  if (dir.exists(adapter_log_dir)) {
    all_files <- list.files(adapter_log_dir)
    cat("  All files in directory:", paste(all_files, collapse=", "), "\n")
  }

  adapter_logs <- list.files(adapter_log_dir, pattern = "^aviti_trimmomatic_pass1_.*\\.log$", full.names = TRUE)
  cat("  Found", length(adapter_logs), "Trimmomatic pass1 stderr logs\n")

  for (log_file in adapter_logs) {
    # Extract sample name from filename (e.g., "aviti_trimmomatic_pass1_KM18-1-AMF_S486.log")
    base_name <- sub("^aviti_trimmomatic_pass1_", "", basename(log_file))
    sample_name <- sub("\\.log$", "", base_name)  # Keep _S## suffix to match DADA2 sample names
    cat("  Parsing:", basename(log_file), "-> sample:", sample_name, "\n")
    stats <- parse_trimmomatic_log(log_file)
    adapter_stats[[sample_name]] <- stats
    cat("    Extracted: input=", stats$input, " output=", stats$output, " pct_removed=", stats$pct_removed, "\n")
  }
} else {
  # Illumina: use cutadapt JSON
  cat("Platform is Illumina: looking for cutadapt JSON files\n")
  if (dir.exists(adapter_log_dir)) {
    all_files <- list.files(adapter_log_dir)
    cat("  All files in directory:", paste(all_files, collapse=", "), "\n")
  }

  adapter_jsons <- list.files(adapter_log_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)
  cat("  Found", length(adapter_jsons), "adapter JSON files matching pattern\n")

  # Also try alternative patterns
  if (length(adapter_jsons) == 0) {
    cat("  Trying alternative pattern *.json...\n")
    adapter_jsons <- list.files(adapter_log_dir, pattern = "\\.json$", full.names = TRUE)
    cat("  Found", length(adapter_jsons), "JSON files with *.json pattern\n")
  }

  for (json_file in adapter_jsons) {
    # Extract sample name, removing cutadapt prefix only (keep _S## suffix to match DADA2 sample names)
    base_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
    sample_name <- base_name  # Keep full name with _S## suffix
    cat("  Parsing:", basename(json_file), "-> sample:", sample_name, "\n")
    stats <- parse_cutadapt_json(json_file)
    adapter_stats[[sample_name]] <- stats
    cat("    Extracted: input=", stats$input, " output=", stats$output, " pct_removed=", stats$pct_removed, "\n")
  }
}

# Parse primer JSONs (cutadapt always for both platforms)
# For both platforms, JSON files are stored in their respective trimmed directories during execution
# but we reference them from .logs/ where the stderr logs are located
if (platform == "aviti") {
  # Aviti primer trimming is done with cutadapt after Trimmomatic adapter pass 1
  # JSON files go to 01b_primer_trimmed/02_logs/ (not deleted before qc runs)
  primer_json_dir <- file.path(base_output_dir, "01b_primer_trimmed", "02_logs")
} else {
  # Illumina primer trimming is done with cutadapt
  # JSON files go to 02_primer_trimmed/02_logs/
  primer_json_dir <- file.path(base_output_dir, "02_primer_trimmed", "02_logs")
}

cat("Looking for primer JSONs in:", primer_json_dir, "\n")
cat("  Directory exists:", dir.exists(primer_json_dir), "\n")

primer_jsons <- list.files(primer_json_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)
cat("  Found", length(primer_jsons), "primer JSON files matching pattern\n")

# Also try alternative patterns
if (length(primer_jsons) == 0) {
  cat("  Trying alternative pattern *.json...\n")
  primer_jsons <- list.files(primer_json_dir, pattern = "\\.json$", full.names = TRUE)
  cat("  Found", length(primer_jsons), "JSON files with *.json pattern\n")
}

primer_stats <- list()
for (json_file in primer_jsons) {
  # Extract sample name, removing cutadapt prefix only (keep _S## suffix to match DADA2 sample names)
  base_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
  sample_name <- base_name  # Keep full name with _S## suffix
  cat("  Parsing:", basename(json_file), "-> sample:", sample_name, "\n")
  stats <- parse_cutadapt_json(json_file)
  primer_stats[[sample_name]] <- stats
  cat("    Extracted: input=", stats$input, " output=", stats$output, " pct_removed=", stats$pct_removed, "\n")
}

cat("\n")
# Debug: show what was parsed
if (length(adapter_stats) > 0) {
  cat("Adapter stats summary:\n")
  for (sname in names(adapter_stats)) {
    cat("  ", sname, ": input=", adapter_stats[[sname]]$input, " output=", adapter_stats[[sname]]$output, " pct_removed=", adapter_stats[[sname]]$pct_removed, "\n")
  }
} else {
  cat("No adapter stats were parsed!\n")
}

if (length(primer_stats) > 0) {
  cat("Primer stats summary:\n")
  for (sname in names(primer_stats)) {
    cat("  ", sname, ": input=", primer_stats[[sname]]$input, " output=", primer_stats[[sname]]$output, " pct_removed=", primer_stats[[sname]]$pct_removed, "\n")
  }
} else {
  cat("No primer stats were parsed!\n")
}

# ==============================================================================
# READ DADA2 SUMMARY
# ==============================================================================
cat("Reading DADA2 summary statistics...\n")

dada2_summary <- read.table(dada2_summary_path, header = TRUE, row.names = 1, sep = "\t")

# ==============================================================================
# COMBINE ALL QC STATS
# ==============================================================================
cat("Combining QC statistics...\n")

# Build combined QC table
sample_names <- rownames(dada2_summary)
qc_df <- data.frame(row.names = sample_names)

# Add reads_start (from adapter stats if available, else from primer stats)
qc_df$reads_start <- NA
for (sname in sample_names) {
  if (sname %in% names(adapter_stats) && !is.na(adapter_stats[[sname]]$input)) {
    qc_df[sname, "reads_start"] <- adapter_stats[[sname]]$input
  } else if (sname %in% names(primer_stats) && !is.na(primer_stats[[sname]]$input)) {
    qc_df[sname, "reads_start"] <- primer_stats[[sname]]$input
  }
}

# Add adapter stats
qc_df$adapter_fwd_pct <- NA
qc_df$adapter_rev_pct <- NA
for (sname in sample_names) {
  if (sname %in% names(adapter_stats)) {
    qc_df[sname, "adapter_fwd_pct"] <- adapter_stats[[sname]]$pct_r1
    qc_df[sname, "adapter_rev_pct"] <- adapter_stats[[sname]]$pct_r2
  }
}


# Add primer stats
qc_df$reads_after_primer <- NA
for (sname in sample_names) {
  if (sname %in% names(primer_stats) && !is.na(primer_stats[[sname]]$output)) {
   qc_df[sname, "reads_after_primer"] <- primer_stats[[sname]]$output
  }
}

qc_df$primer_fwd_pct <- NA
qc_df$primer_rev_pct <- NA
qc_df$primer_pct_removed <- NA
for (sname in sample_names) {
  if (sname %in% names(primer_stats)) {
    qc_df[sname, "primer_fwd_pct"] <- primer_stats[[sname]]$pct_r1
    qc_df[sname, "primer_rev_pct"] <- primer_stats[[sname]]$pct_r2
    qc_df[sname, "primer_pct_removed"] <- primer_stats[[sname]]$pct_removed
  }
}

# Add DADA2 stats
qc_df$reads_dada2_input <- dada2_summary[sample_names, "dada2_input"]

qc_df$dada2_filtered <- dada2_summary[sample_names, "filtered"]
qc_df$pct_filtered_out <- ifelse(qc_df$reads_dada2_input > 0,
  100 * (qc_df$reads_dada2_input - qc_df$dada2_filtered) / qc_df$reads_dada2_input, NA)

qc_df$dada_f <- dada2_summary[sample_names, "dada_f"]
qc_df$dada_r <- dada2_summary[sample_names, "dada_r"]
qc_df$dada2_merged <- dada2_summary[sample_names, "merged"]
qc_df$merge_pct <- ifelse(qc_df$dada2_filtered > 0, 100 * qc_df$dada2_merged / qc_df$dada2_filtered, NA)

qc_df$dada2_nonchim <- dada2_summary[sample_names, "nonchim"]
qc_df$pct_chimeric <- ifelse(qc_df$dada2_merged > 0,
  100 * (qc_df$dada2_merged - qc_df$dada2_nonchim) / qc_df$dada2_merged, NA)


# Placeholder for taxonomy assignment percentage (will be filled in from taxonomy file if available)
qc_df$taxonomy_assigned_pct <- NA

# Find and read taxonomy file to calculate taxonomy_assigned_pct
taxa_files <- list.files(dada2_dir, pattern = "__combined_sequences_ASVtaxa_.*\\.txt$", full.names = TRUE)
taxa_files <- taxa_files[!grepl("bootstrap", taxa_files)]
if (length(taxa_files) == 1) {
  tryCatch({
    taxa_table <- read.table(taxa_files[1], header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
    first_tax_col <- intersect(c("Kingdom", "Domain"), colnames(taxa_table))[1]
    if (!is.na(first_tax_col)) {
      asv_tax_pct <- 100 * sum(!is.na(taxa_table[[first_tax_col]])) / nrow(taxa_table)
      # Set same value for all samples (per-run metric)
      qc_df$taxonomy_assigned_pct <- asv_tax_pct
    }
  }, error = function(e) {
    cat("Warning: Could not calculate taxonomy_assigned_pct:", e$message, "\n")
  })
}

# ==============================================================================
# WRITE QC SUMMARY FILE
# ==============================================================================
cat("Writing QC summary file...\n")

qc_summary_path <- file.path(qc_outdir, "qc_summary.txt")

# Open file for writing header comments
sink(qc_summary_path)

cat("# QC Summary - DADA2 Pipeline\n")
cat("# Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("# Amplicon: ", amplicon, "\n")
cat("#\n")
cat("# COLUMN DESCRIPTIONS:\n")
cat("#   sample:                   Sample name (or 'average' for row means)\n")
cat("#   reads_start:              Total input reads\n")
cat("#   adapter_fwd_pct:          % of forward reads with adapter sequence detected\n")
cat("#   adapter_rev_pct:          % of reverse reads with adapter sequence detected\n")
cat("#   reads_after_primer:       Reads after primer trimming\n")
cat("#   primer_fwd_pct:           % of forward reads with primer detected\n")
cat("#   primer_rev_pct:           % of reverse reads with primer detected\n")
cat("#   primer_pct_removed:       % reads removed by primer trimming (--discard-untrimmed)\n")
cat("#   reads_dada2_input:        Reads entering DADA2 denoising\n")
cat("#   dada2_filtered:           Reads after DADA2 quality filtering\n")
cat("#   pct_filtered_out:         % reads removed by DADA2 quality filtering\n")
cat("#   dada_f:                   total forward reads\n")
cat("#   dada_r:                   total reverse reads\n")
cat("#   dada2_merged:             number of merged reads\n")
cat("#   merge_pct:                % of filtered reads that merged successfully\n")
cat("#   dada2_nonchim:            number of non-chimeric reads\n")
cat("#   pct_chimeric:             % of merged reads that were chimeric\n")
cat("#   taxonomy_assigned_pct:    % of ASVs with taxonomy assignment\n")
cat("#\n")
cat("# EXPECTED VALUES (approximate):\n")
cat("#   adapter_fwd/rev_pct:      variable (low for long amplicons >300bp; high for 16S/short)\n")
cat("#   primer_fwd/rev_pct:       typically 70-95% detected (18S-V4 R2 often 50-70%)\n")
cat("#   primer_pct_removed:       typically 5-25% removed\n")
cat("#   pct_filtered_out:         typically 15-30% removed (70-85% retained)\n")
cat("#   merge_pct:                16S: 75-95%, ITS: 60-85%, 18S-V4: 30-70%\n")
cat("#   pct_chimeric:             typically 2-15% (85-98% retained)\n")
cat("#   taxonomy_assigned_pct:    typically >90% if database matches samples\n")
cat("#\n")

sink()

# Convert rownames to a "sample" column
qc_df$sample <- rownames(qc_df)
qc_df <- qc_df[, c("sample", setdiff(colnames(qc_df), "sample"))]  # Move sample to first column

# Calculate averages for all numeric columns
avg_row <- data.frame(sample = "average")
for (col in colnames(qc_df)[-1]) {  # Skip "sample" column
  if (is.numeric(qc_df[[col]])) {
    avg_row[[col]] <- mean(qc_df[[col]], na.rm = TRUE)
  } else {
    avg_row[[col]] <- NA
  }
}

# Ensure column order matches qc_df
avg_row <- avg_row[, colnames(qc_df)]

# Create blank row with proper column types (copy structure from qc_df, then set all to NA)
blank_row <- qc_df[1, , drop = FALSE]
blank_row[,] <- NA

# Add blank row and average row
qc_df_with_avg <- rbind(qc_df, blank_row, stringsAsFactors = FALSE)
qc_df_with_avg <- rbind(qc_df_with_avg, avg_row, stringsAsFactors = FALSE)

# Write the table with sample as a regular column (na="" makes blank row truly blank)
write.table(qc_df_with_avg,
  file = qc_summary_path,
  sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE, row.names = FALSE, na = ""
)

cat("QC summary written to:", qc_summary_path, "\n\n")

# ==============================================================================
# GENERATE FIGURES
# ==============================================================================
cat("Generating QC figures...\n\n")

# Quality profiles
cat("Generating quality profile figures...\n")

# Read filtered reads for quality profile
primer_fqs_fwd <- sort(list.files(primer_trimmed_dir, pattern = "_R1_001.fastq.gz", full.names = TRUE))
primer_fqs_rev <- sort(list.files(primer_trimmed_dir, pattern = "_R2_001.fastq.gz", full.names = TRUE))

cat("  Found ", length(primer_fqs_fwd), " forward and ", length(primer_fqs_rev), " reverse fastq files\n")

# Limit to first 3 samples for plotting
n_samples <- min(3, length(primer_fqs_fwd))
if (n_samples > 0) {
  cat("  Processing quality profiles for first", n_samples, "samples...\n")

  tryCatch({
    cat("    Plotting forward reads quality profile...\n")
    p <- plotQualityProfile(primer_fqs_fwd[1:n_samples])
    ggsave(file.path(qc_outdir, "figures", "quality_profile_fwd.pdf"), plot = p, width = 12, height = 7)
    cat("  ✓ quality_profile_fwd.pdf\n")
  }, error = function(e) {
    cat("  ✗ quality_profile_fwd.pdf:", e$message, "\n")
  })

  tryCatch({
    cat("    Plotting reverse reads quality profile...\n")
    p <- plotQualityProfile(primer_fqs_rev[1:n_samples])
    ggsave(file.path(qc_outdir, "figures", "quality_profile_rev.pdf"), plot = p, width = 12, height = 7)
    cat("  ✓ quality_profile_rev.pdf\n")
  }, error = function(e) {
    cat("  ✗ quality_profile_rev.pdf:", e$message, "\n")
  })
} else {
  cat("  Warning: No fastq files found in", primer_trimmed_dir, "\n")
}

# Error model plots
cat("Generating error model figures...\n")

tryCatch({
  cat("  Loading error models...\n")
  errF_path <- file.path(dada2_dir, "errF.rds")
  errR_path <- file.path(dada2_dir, "errR.rds")

  if (!file.exists(errF_path)) {
    cat("  Warning: errF.rds not found at", errF_path, "\n")
  } else if (!file.exists(errR_path)) {
    cat("  Warning: errR.rds not found at", errR_path, "\n")
  } else {
    errF <- readRDS(errF_path)
    errR <- readRDS(errR_path)

    tryCatch({
      cat("    Plotting forward error model...\n")
      p <- plotErrors(errF, nominalQ = TRUE)
      ggsave(file.path(qc_outdir, "figures", "error_model_fwd.pdf"), plot = p, width = 12, height = 7)
      cat("  ✓ error_model_fwd.pdf\n")
    }, error = function(e) {
      cat("  ✗ error_model_fwd.pdf:", e$message, "\n")
    })

    # Only plot reverse error model if errR is not NULL (i.e., not in forward-only mode)
    if (!is.null(errR)) {
      tryCatch({
        cat("    Plotting reverse error model...\n")
        p <- plotErrors(errR, nominalQ = TRUE)
        ggsave(file.path(qc_outdir, "figures", "error_model_rev.pdf"), plot = p, width = 12, height = 7)
        cat("  ✓ error_model_rev.pdf\n")
      }, error = function(e) {
        cat("  ✗ error_model_rev.pdf:", e$message, "\n")
      })
    } else {
      cat("    Reverse error model skipped (forward-only mode)\n")
    }
  }
}, error = function(e) {
  cat("  Error loading error models:", e$message, "\n")
})

# Sequence length distribution
cat("Generating sequence length distribution figure...\n")

tryCatch({
  seqtab_path <- file.path(dada2_dir, "seqtab_nochim.rds")
  if (!file.exists(seqtab_path)) {
    cat("  Warning: seqtab_nochim.rds not found at", seqtab_path, "\n")
  } else {
    seqtab.nochim <- readRDS(seqtab_path)
    seq_lengths <- nchar(colnames(seqtab.nochim))
    cat("  ", length(seq_lengths), "ASVs found, length range:", min(seq_lengths), "-", max(seq_lengths), "bp\n")

    tryCatch({
      pdf(file.path(qc_outdir, "figures", "sequence_length_distribution.pdf"), width = 10, height = 6)
      hist(seq_lengths, breaks = 50, main = "Sequence Length Distribution", xlab = "Length (bp)", ylab = "Frequency")
      dev.off()
      cat("  ✓ sequence_length_distribution.pdf\n")
    }, error = function(e) {
      cat("  ✗ sequence_length_distribution.pdf:", e$message, "\n")
    })
  }
}, error = function(e) {
  cat("  Error processing sequence table:", e$message, "\n")
})

# Read tracking
cat("Generating read tracking figure...\n")

tryCatch({
  # Define all possible steps in order
  step_defs <- list(
    list(name = "Start",          col = "reads_start"),
    list(name = "After Primer",   col = "reads_after_primer"),
    list(name = "DADA2 Input",    col = "reads_dada2_input"),
    list(name = "After Filter",   col = "dada2_filtered"),
    list(name = "After Merge",    col = "dada2_merged"),
    list(name = "After Chimera",  col = "dada2_nonchim")
  )

  # Build step_data with parallel step_names_used vector
  step_data <- data.frame(step_name = character(), mean_reads = numeric())
  step_names_used <- character()

  for (s in step_defs) {
    if (s$col %in% colnames(qc_df) && any(!is.na(qc_df[[s$col]]))) {
      vals <- qc_df[[s$col]]
      step_data <- rbind(step_data, data.frame(step_name = s$name, mean_reads = mean(vals, na.rm = TRUE)))
      step_names_used <- c(step_names_used, s$name)
    }
  }

  if (nrow(step_data) > 0) {
    # Create ggplot with better spacing, using step_names_used for proper factor levels
    p <- ggplot(step_data, aes(x = factor(step_name, levels = step_names_used), y = mean_reads, fill = step_name)) +
      geom_bar(stat = "identity", color = "black", linewidth = 0.7) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11, margin = margin(t = 10)),
        axis.text.y = element_text(size = 11, margin = margin(r = 10)),
        axis.title = element_text(size = 12, margin = margin(t = 10, r = 10)),
        plot.title = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)),
        legend.position = "none"
      ) +
      labs(
        title = "Mean Reads per Sample at Each Pipeline Step",
        x = "Pipeline Step",
        y = "Mean Number of Reads"
      ) +
      scale_fill_brewer(palette = "Set3")

    ggsave(file.path(qc_outdir, "figures", "read_tracking.pdf"), plot = p, width = 12, height = 7)
    cat("  read_tracking.pdf\n")
  } else {
    cat("  Warning: could not generate read_tracking.pdf - no data available\n")
  }
}, error = function(e) {
  cat("  Warning: could not generate read_tracking.pdf:", e$message, "\n")
})

cat("\nQC figures generated in:", file.path(qc_outdir, "figures"), "\n")

# ==============================================================================
# GENERATE README FOR FIGURES
# ==============================================================================
readme_path <- file.path(qc_outdir, "figures", "README.md")
sink(readme_path)

cat("# QC Figures - DADA2 Pipeline\n\n")
cat("This directory contains quality control visualizations for the DADA2 amplicon processing pipeline.\n\n")

cat("## Files\n\n")

cat("### quality_profile_fwd.pdf\n")
cat("**Quality scores across forward (R1) reads**\n\n")
cat("Shows the PHRED quality score distribution for the first 3 samples (up to 100,000 reads per sample):\n")
cat("- **X-axis**: Position along the read (base pair position)\n")
cat("- **Y-axis**: Quality score (Q-score, range 0-40)\n")
cat("- **Gray shading**: Shows the quartiles (25th, 50th, 75th percentiles) of quality at each position\n")
cat("- **Red line**: Mean quality score\n\n")
cat("**Interpretation**:\n")
cat("- Quality should remain **>30** (1 error per 1000 bases) through most of the read\n")
cat("- Quality drops toward the read end (normal)\n")
cat("- If quality drops sharply mid-read, suggests sequencing issues\n")
cat("- For paired-end reads, forward reads are typically higher quality than reverse\n\n")

cat("### quality_profile_rev.pdf\n")
cat("**Quality scores across reverse (R2) reads**\n\n")
cat("Same interpretation as quality_profile_fwd.pdf, but for reverse reads.\n\n")
cat("**Interpretation**:\n")
cat("- Reverse reads typically show lower quality than forward reads\n")
cat("- This is expected and normal for Illumina sequencing\n")
cat("- DADA2 quality filtering will trim or discard low-quality reads\n\n")

cat("### error_model_fwd.pdf\n")
cat("**Error model for forward reads**\n\n")
cat("Shows the relationship between observed quality score (Q) and actual error rate:\n")
cat("- **X-axis**: Nominal quality score (Q-score)\n")
cat("- **Y-axis**: Observed error rate (log scale)\n")
cat("- **Black dots**: Observed errors at each quality score\n")
cat("- **Red line**: Estimated error rate (loess fit)\n")
cat("- **Blue line**: Expected error rate (ideal model)\n\n")
cat("**Interpretation**:\n")
cat("- The red line should closely follow the blue line\n")
cat("- If red line is far above blue line at high Q-scores, suggests systematic errors\n")
cat("- Deviations at low Q-scores are less concerning\n")
cat("- This error model is used by DADA2 for denoising\n\n")

cat("### error_model_rev.pdf\n")
cat("**Error model for reverse reads**\n\n")
cat("Same interpretation as error_model_fwd.pdf, but for reverse reads.\n\n")

cat("### sequence_length_distribution.pdf\n")
cat("**Distribution of final ASV (Amplicon Sequence Variant) lengths**\n\n")
cat("Histogram showing the length (in base pairs) of all unique sequences in the final ASV table:\n")
cat("- **X-axis**: Sequence length (bp)\n")
cat("- **Y-axis**: Number of ASVs\n")
cat("- **Title shows**: Total ASVs, length range (min-max)\n\n")
cat("**Interpretation**:\n")
cat("- Should show a narrow, tight distribution around the expected amplicon length\n")
cat("- For 16S-V4: expect ~400 bp\n")
cat("- For ITS: expect wider distribution (200-600 bp)\n")
cat("- For 18S-AMF: expect ~300-550 bp\n")
cat("- Outliers (very short or very long) suggest:\n")
cat("  - Short: primer dimers or non-specific amplification\n")
cat("  - Long: chimeras or multiple amplicons\n")
cat("- Multiple peaks suggest mixed amplicon types or contamination\n\n")

cat("### read_tracking.pdf\n")
cat("**Read counts at each pipeline step**\n\n")
cat("Bar chart showing the mean number of reads retained at each processing step:\n")
cat("- **X-axis**: Pipeline steps (Start → Adapter → Primer → DADA2 Input → Filter → Merge → Chimera)\n")
cat("- **Y-axis**: Mean number of reads (log scale may be used if range is very large)\n")
cat("- **Bar color**: Different color for each step\n\n")
cat("**Interpretation** (expected % retention between steps):\n")
cat("- **Adapter trimming**: 98-100% (trimming doesn't remove reads for long amplicons)\n")
cat("- **Primer trimming**: 85-97% (removes reads without both primers)\n")
cat("- **Quality filtering**: 70-85% (removes low-quality sequences)\n")
cat("- **Merging**: 30-95% depending on amplicon (requires sufficient overlap)\n")
cat("- **Chimera removal**: 85-98% (removes chimeric sequences)\n\n")
cat("**Red flags**:\n")
cat("- Large drop at primer step: wrong primers or strand orientation\n")
cat("- Minimal merge rate: primers too far apart or low overlap\n")
cat("- Unexpectedly low values at any step: check raw data quality\n\n")

cat("## Summary Statistics\n\n")
cat("For detailed QC statistics by sample, see `../qc_summary.txt`.\n\n")
cat("Key metrics to monitor:\n")
cat("- **pct_filtered_out**: % reads removed by quality filtering (15-30% expected)\n")
cat("- **merge_pct**: % of filtered reads that merged (30-95% depending on amplicon)\n")
cat("- **pct_chimeric**: % of merged reads that were chimeric (2-15% expected)\n")
cat("- **taxonomy_assigned_pct**: % of ASVs with successful taxonomy assignment (>90% expected)\n\n")

cat("## Troubleshooting\n\n")
cat("**Quality profiles show consistently low quality?**\n")
cat("- Check your input fastq files for sequencing issues\n")
cat("- Consider using more stringent quality filtering in Snakefile\n\n")

cat("**Error model doesn't fit well?**\n")
cat("- Can happen with unusual sample types or very high error rates\n")
cat("- DADA2 will still work, but accuracy may be lower\n\n")

cat("**Sequence length distribution shows multiple peaks?**\n")
cat("- Check for primer dimers or non-specific amplification\n")
cat("- Review primer sequences and alignment\n\n")

cat("**Low merge rate?**\n")
cat("- Primers may be too far apart for overlap\n")
cat("- Read quality may be too low toward the ends\n")
cat("- For 16S, 75-95% merge rate is typical; for ITS, 60-85%\n\n")

cat("**High chimera rate (>15%)?**\n")
cat("- Suggests overamplification or non-specific primer binding\n")
cat("- Check PCR conditions and primer specificity\n")

sink()

cat("QC README generated at:", readme_path, "\n")
cat("\nQC figures generated in:", file.path(qc_outdir, "figures"), "\n")
cat("\nDone.\n")
