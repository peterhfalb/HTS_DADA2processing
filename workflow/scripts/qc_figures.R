#!/usr/bin/env Rscript
# QC figures and stats generation for DADA2 pipeline
# Usage: Rscript qc_figures.R <qc_outdir> <amplicon> <dada2_summary_path> <dada2_dir> <primer_trimmed_dir>

library(dada2)

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
# PARSE CUTADAPT LOGS
# ==============================================================================
cat("Parsing cutadapt logs...\n")

parse_cutadapt_log <- function(log_path) {
  tryCatch({
    lines <- readLines(log_path)

    # Extract total read pairs processed
    total_line <- lines[grepl("Total read pairs processed", lines)]
    if (length(total_line) == 0) return(list(input = NA, passing = NA))

    total <- as.integer(gsub("[^0-9]", "", strsplit(total_line, ":", fixed = TRUE)[[1]][2]))

    # Extract passing read pairs
    passing_line <- lines[grepl("Pairs written.*passing", lines)]
    if (length(passing_line) == 0) return(list(input = total, passing = NA))

    passing_str <- strsplit(passing_line, ":", fixed = TRUE)[[1]][2]
    passing <- as.integer(gsub("[^0-9]", "", strsplit(passing_str, " ", fixed = TRUE)[[1]][1]))

    list(input = total, passing = passing)
  }, error = function(e) {
    cat("Warning: could not parse log", log_path, "\n")
    list(input = NA, passing = NA)
  })
}

# Parse adapter logs
adapter_log_dir <- file.path(dirname(dirname(primer_trimmed_dir)), "01_adapter", "01_logs")
adapter_logs <- list.files(adapter_log_dir, pattern = "^cutadapt\\..*\\.log\\.txt$", full.names = TRUE)

adapter_stats <- list()
for (log_file in adapter_logs) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.log\\.txt$", "", basename(log_file)))
  stats <- parse_cutadapt_log(log_file)
  adapter_stats[[sample_name]] <- stats
}

# Parse primer logs
primer_log_dir <- file.path(primer_trimmed_dir, "02_logs")
primer_logs <- list.files(primer_log_dir, pattern = "^cutadapt\\..*\\.log\\.txt$", full.names = TRUE)

primer_stats <- list()
for (log_file in primer_logs) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.log\\.txt$", "", basename(log_file)))
  stats <- parse_cutadapt_log(log_file)
  primer_stats[[sample_name]] <- stats
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

# Add adapter stats
qc_df$adapter_input <- NA
qc_df$adapter_passing <- NA
for (sname in sample_names) {
  if (sname %in% names(adapter_stats)) {
    qc_df[sname, "adapter_input"] <- adapter_stats[[sname]]$input
    qc_df[sname, "adapter_passing"] <- adapter_stats[[sname]]$passing
  }
}

# Add primer stats
qc_df$primer_input <- NA
qc_df$primer_passing <- NA
for (sname in sample_names) {
  if (sname %in% names(primer_stats)) {
    qc_df[sname, "primer_input"] <- primer_stats[[sname]]$input
    qc_df[sname, "primer_passing"] <- primer_stats[[sname]]$passing
  }
}

# Add DADA2 stats
qc_df$dada2_input <- dada2_summary[sample_names, "dada2_input"]
qc_df$filtered <- dada2_summary[sample_names, "filtered"]
qc_df$dada_f <- dada2_summary[sample_names, "dada_f"]
qc_df$dada_r <- dada2_summary[sample_names, "dada_r"]
qc_df$merged <- dada2_summary[sample_names, "merged"]
qc_df$nonchim <- dada2_summary[sample_names, "nonchim"]

# Calculate percentages
qc_df$adapter_pct <- 100 * qc_df$adapter_passing / qc_df$adapter_input
qc_df$primer_pct <- 100 * qc_df$primer_passing / qc_df$primer_input
qc_df$dada2_filter_pct <- 100 * qc_df$filtered / qc_df$dada2_input
qc_df$merge_pct <- 100 * qc_df$merged / qc_df$filtered
qc_df$chimera_pct <- 100 * qc_df$nonchim / qc_df$merged

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
cat("# EXPECTED VALUES (approximate):\n")
cat("#   adapter_retention:    typically 85-98% (low = adapter contamination or poor quality)\n")
cat("#   primer_retention:     typically 70-95% (low = wrong primer or orientation issues)\n")
cat("#   dada2_filter:         typically 70-85% (< 50% suggests quality problems)\n")
cat("#   merge_rate:           typically 75-95% (low = insufficient overlap or high error rate)\n")
cat("#   chimera_removal:      typically 85-98% remaining (> 15% chimeric = concern)\n")
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

# Limit to first 3 samples for plotting
n_samples <- min(3, length(primer_fqs_fwd))
if (n_samples > 0) {
  tryCatch({
    pdf(file.path(qc_outdir, "figures", "quality_profile_fwd.pdf"), width = 10, height = 6)
    plotQualityProfile(primer_fqs_fwd[1:n_samples])
    dev.off()
    cat("  quality_profile_fwd.pdf\n")
  }, error = function(e) {
    cat("  Warning: could not generate quality_profile_fwd.pdf:", e$message, "\n")
  })

  tryCatch({
    pdf(file.path(qc_outdir, "figures", "quality_profile_rev.pdf"), width = 10, height = 6)
    plotQualityProfile(primer_fqs_rev[1:n_samples])
    dev.off()
    cat("  quality_profile_rev.pdf\n")
  }, error = function(e) {
    cat("  Warning: could not generate quality_profile_rev.pdf:", e$message, "\n")
  })
}

# Error model plots
cat("Generating error model figures...\n")

errF <- readRDS(file.path(dada2_dir, "errF.rds"))
errR <- readRDS(file.path(dada2_dir, "errR.rds"))

tryCatch({
  pdf(file.path(qc_outdir, "figures", "error_model_fwd.pdf"), width = 10, height = 6)
  plotErrors(errF, nominalQ = TRUE)
  dev.off()
  cat("  error_model_fwd.pdf\n")
}, error = function(e) {
  cat("  Warning: could not generate error_model_fwd.pdf:", e$message, "\n")
})

tryCatch({
  pdf(file.path(qc_outdir, "figures", "error_model_rev.pdf"), width = 10, height = 6)
  plotErrors(errR, nominalQ = TRUE)
  dev.off()
  cat("  error_model_rev.pdf\n")
}, error = function(e) {
  cat("  Warning: could not generate error_model_rev.pdf:", e$message, "\n")
})

# Sequence length distribution
cat("Generating sequence length distribution figure...\n")

seqtab.nochim <- readRDS(file.path(dada2_dir, "seqtab_nochim.rds"))
seq_lengths <- nchar(colnames(seqtab.nochim))

tryCatch({
  pdf(file.path(qc_outdir, "figures", "sequence_length_distribution.pdf"), width = 10, height = 6)
  hist(seq_lengths, breaks = 50, main = "Sequence Length Distribution", xlab = "Length (bp)", ylab = "Frequency")
  dev.off()
  cat("  sequence_length_distribution.pdf\n")
}, error = function(e) {
  cat("  Warning: could not generate sequence_length_distribution.pdf:", e$message, "\n")
})

# Read tracking
cat("Generating read tracking figure...\n")

tryCatch({
  # Calculate mean reads at each step
  steps <- c("adapter_input", "adapter_passing", "primer_passing", "dada2_input", "filtered", "merged", "nonchim")
  step_data <- matrix(NA, nrow = length(sample_names), ncol = length(steps))

  for (i in seq_along(steps)) {
    if (steps[i] %in% colnames(qc_df)) {
      step_data[, i] <- qc_df[sample_names, steps[i]]
    }
  }

  # Use mean reads for visualization
  mean_reads <- colMeans(step_data, na.rm = TRUE)

  pdf(file.path(qc_outdir, "figures", "read_tracking.pdf"), width = 12, height = 6)

  # Create barplot with proper colors
  colors <- rainbow(length(steps))
  barplot(mean_reads, names.arg = steps, main = "Mean Reads per Sample at Each Pipeline Step",
    ylab = "Mean Number of Reads", xlab = "Pipeline Step", col = colors, las = 2, cex.names = 0.8)

  dev.off()
  cat("  read_tracking.pdf\n")
}, error = function(e) {
  cat("  Warning: could not generate read_tracking.pdf:", e$message, "\n")
})

cat("\nQC figures generated in:", file.path(qc_outdir, "figures"), "\n")
cat("\nDone.\n")
