#!/usr/bin/env Rscript
# QC figures and stats generation for DADA2 pipeline
# Usage: Rscript qc_figures.R <qc_outdir> <amplicon> <dada2_summary_path> <dada2_dir> <primer_trimmed_dir>

library(dada2)
library(jsonlite)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript qc_figures.R <qc_outdir> <amplicon> <dada2_summary_path> <dada2_dir> <primer_trimmed_dir>")
}

qc_outdir <- args[1]
amplicon <- args[2]
dada2_summary_path <- args[3]
dada2_dir <- args[4]
primer_trimmed_dir <- args[5]

cat("QC Figures Generation\n")
cat("====================\n")
cat("QC output dir:        ", qc_outdir, "\n")
cat("Amplicon type:        ", amplicon, "\n")
cat("DADA2 summary:        ", dada2_summary_path, "\n")
cat("DADA2 dir:            ", dada2_dir, "\n")
cat("Primer-trimmed dir:   ", primer_trimmed_dir, "\n\n")

# Ensure output directories exist
dir.create(file.path(qc_outdir, "figures"), showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# PARSE CUTADAPT JSON FILES
# ==============================================================================
cat("Parsing cutadapt JSON files...\n")

parse_cutadapt_json <- function(json_path) {
  tryCatch({
    data <- jsonlite::fromJSON(json_path)

    # Extract read counts - handle different cutadapt versions
    if (!is.null(data$read_counts)) {
      rc <- data$read_counts
    } else if (!is.null(data$statistics$input)) {
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
    cat("Warning: could not parse JSON", json_path, ":", e$message, "\n")
    list(input = NA, output = NA, r1_with_adapter = NA, r2_with_adapter = NA, pct_r1 = NA, pct_r2 = NA, pct_removed = NA)
  })
}

# Parse adapter JSONs
adapter_json_dir <- file.path(dirname(dirname(primer_trimmed_dir)), "01_adapter", "01_logs")
cat("Looking for adapter JSONs in:", adapter_json_dir, "\n")
adapter_jsons <- list.files(adapter_json_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)
cat("  Found", length(adapter_jsons), "adapter JSON files\n")

adapter_stats <- list()
for (json_file in adapter_jsons) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
  stats <- parse_cutadapt_json(json_file)
  adapter_stats[[sample_name]] <- stats
}

# Parse primer JSONs
primer_json_dir <- file.path(primer_trimmed_dir, "02_logs")
cat("Looking for primer JSONs in:", primer_json_dir, "\n")
primer_jsons <- list.files(primer_json_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)
cat("  Found", length(primer_jsons), "primer JSON files\n")

primer_stats <- list()
for (json_file in primer_jsons) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
  stats <- parse_cutadapt_json(json_file)
  primer_stats[[sample_name]] <- stats
}

# Debug: show what was parsed
if (length(adapter_stats) > 0) {
  cat("Adapter stats summary:\n")
  for (sname in names(adapter_stats)[1:min(2, length(adapter_stats))]) {
    cat("  ", sname, ": input=", adapter_stats[[sname]]$input, " output=", adapter_stats[[sname]]$output, " pct_removed=", adapter_stats[[sname]]$pct_removed, "\n")
  }
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
qc_df$adapter_pct_removed <- NA
for (sname in sample_names) {
  if (sname %in% names(adapter_stats)) {
    qc_df[sname, "adapter_fwd_pct"] <- adapter_stats[[sname]]$pct_r1
    qc_df[sname, "adapter_rev_pct"] <- adapter_stats[[sname]]$pct_r2
    qc_df[sname, "adapter_pct_removed"] <- adapter_stats[[sname]]$pct_removed
  }
}

# Add primer stats
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
qc_df$dada_f <- dada2_summary[sample_names, "dada_f"]
qc_df$dada_r <- dada2_summary[sample_names, "dada_r"]
qc_df$dada2_merged <- dada2_summary[sample_names, "merged"]
qc_df$dada2_nonchim <- dada2_summary[sample_names, "nonchim"]

# Calculate percentage metrics
qc_df$dada2_filter_pct <- ifelse(qc_df$reads_dada2_input > 0, 100 * qc_df$dada2_filtered / qc_df$reads_dada2_input, NA)
qc_df$derep_pct_removed <- NA  # Placeholder - would need dereplicated counts from DADA2
qc_df$denoise_pct_removed <- ifelse(qc_df$dada_f > 0 | qc_df$dada_r > 0,
                                     100 * (qc_df$dada2_filtered - pmax(qc_df$dada_f, qc_df$dada_r, na.rm=TRUE)) / qc_df$dada2_filtered, NA)
qc_df$merge_pct <- ifelse(qc_df$dada2_filtered > 0, 100 * qc_df$dada2_merged / qc_df$dada2_filtered, NA)
qc_df$chimera_pct <- ifelse(qc_df$dada2_merged > 0, 100 * qc_df$dada2_nonchim / qc_df$dada2_merged, NA)

# Placeholder for taxonomy assignment percentage (would need actual taxonomy data)
qc_df$taxonomy_assigned_pct <- NA

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
cat("#   reads_start:              Total input reads\n")
cat("#   adapter_fwd_pct:          % of forward reads with adapter detected\n")
cat("#   adapter_rev_pct:          % of reverse reads with adapter detected\n")
cat("#   adapter_pct_removed:      % reads removed by adapter trimming\n")
cat("#   primer_fwd_pct:           % of forward reads with primer detected\n")
cat("#   primer_rev_pct:           % of reverse reads with primer detected\n")
cat("#   primer_pct_removed:       % reads removed by primer trimming\n")
cat("#   reads_dada2_input:        Reads entering DADA2 denoising\n")
cat("#   dada2_filter_pct:         % reads retained by quality filtering\n")
cat("#   denoise_pct_removed:      % reads removed by denoising\n")
cat("#   merge_pct:                % of denoised reads that merged successfully\n")
cat("#   chimera_pct:              % of merged reads retained after chimera removal\n")
cat("#   taxonomy_assigned_pct:    % of ASVs with taxonomy assignment\n")
cat("#\n")
cat("# EXPECTED VALUES (approximate):\n")
cat("#   adapter_fwd/rev_pct:      typically 85-98% detected\n")
cat("#   adapter_pct_removed:      typically 2-15% removed\n")
cat("#   primer_fwd/rev_pct:       typically 70-95% detected (18S-V4 R2 often 50-70%)\n")
cat("#   primer_pct_removed:       typically 5-25% removed\n")
cat("#   dada2_filter_pct:         typically 70-85% retained (< 50% suggests quality issues)\n")
cat("#   merge_pct:                16S: 75-95%, ITS: 60-85%, 18S-V4: 30-70%\n")
cat("#   chimera_pct:              typically 85-98% retained (> 15% chimeric = concern)\n")
cat("#\n")

sink()

# Write the table
write.table(qc_df,
  file = qc_summary_path,
  sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE, row.names = TRUE
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

    tryCatch({
      cat("    Plotting reverse error model...\n")
      p <- plotErrors(errR, nominalQ = TRUE)
      ggsave(file.path(qc_outdir, "figures", "error_model_rev.pdf"), plot = p, width = 12, height = 7)
      cat("  ✓ error_model_rev.pdf\n")
    }, error = function(e) {
      cat("  ✗ error_model_rev.pdf:", e$message, "\n")
    })
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
  library(ggplot2)

  # Build read tracking data from all available columns
  step_labels <- c("reads_start", "adapter_passing", "primer_passing", "reads_dada2_input", "dada2_filtered", "dada2_merged", "dada2_nonchim")
  step_names <- c("Start", "Adapter\nRemoved", "Primer\nRemoved", "DADA2\nInput", "After\nFilter", "After\nMerge", "After\nChimera")

  step_data <- data.frame(step_label = character(), step_name = character(), mean_reads = numeric())

  for (i in seq_along(step_labels)) {
    col <- step_labels[i]

    # Special handling for adapter_passing and primer_passing (calculate from input - removed)
    if (col == "adapter_passing" && "reads_start" %in% colnames(qc_df) && "adapter_pct_removed" %in% colnames(qc_df)) {
      vals <- qc_df$reads_start * (100 - qc_df$adapter_pct_removed) / 100
    } else if (col == "primer_passing" && "adapter_passing" %in% colnames(qc_df) && "primer_pct_removed" %in% colnames(qc_df)) {
      # Use previously calculated adapter_passing values
      adapter_vals <- qc_df$reads_start * (100 - qc_df$adapter_pct_removed) / 100
      vals <- adapter_vals * (100 - qc_df$primer_pct_removed) / 100
    } else if (col %in% colnames(qc_df)) {
      vals <- qc_df[[col]]
    } else {
      vals <- NA
    }

    if (any(!is.na(vals))) {
      means <- mean(vals, na.rm = TRUE)
      step_data <- rbind(step_data, data.frame(step_label = col, step_name = step_names[i], mean_reads = means))
    }
  }

  if (nrow(step_data) > 0) {
    # Create ggplot with better spacing
    p <- ggplot(step_data, aes(x = factor(step_name, levels = step_names[1:nrow(step_data)]), y = mean_reads, fill = step_name)) +
      geom_bar(stat = "identity", color = "black", size = 0.7) +
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
cat("\nDone.\n")
