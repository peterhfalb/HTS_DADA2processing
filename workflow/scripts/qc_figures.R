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

    # Extract read counts
    total_reads <- data$read_counts$input
    reads_with_adapter_r1 <- data$read_counts$read1_with_adapter
    reads_with_adapter_r2 <- data$read_counts$read2_with_adapter
    output_reads <- data$read_counts$output

    list(
      input = total_reads,
      output = output_reads,
      r1_with_adapter = reads_with_adapter_r1,
      r2_with_adapter = reads_with_adapter_r2,
      pct_r1 = ifelse(!is.null(reads_with_adapter_r1), 100 * reads_with_adapter_r1 / total_reads, NA),
      pct_r2 = ifelse(!is.null(reads_with_adapter_r2), 100 * reads_with_adapter_r2 / total_reads, NA)
    )
  }, error = function(e) {
    cat("Warning: could not parse JSON", json_path, "\n")
    list(input = NA, output = NA, r1_with_adapter = NA, r2_with_adapter = NA, pct_r1 = NA, pct_r2 = NA)
  })
}

# Parse adapter JSONs
adapter_json_dir <- file.path(dirname(dirname(primer_trimmed_dir)), "01_adapter", "01_logs")
adapter_jsons <- list.files(adapter_json_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)

adapter_stats <- list()
for (json_file in adapter_jsons) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
  stats <- parse_cutadapt_json(json_file)
  adapter_stats[[sample_name]] <- stats
}

# Parse primer JSONs
primer_json_dir <- file.path(primer_trimmed_dir, "02_logs")
primer_jsons <- list.files(primer_json_dir, pattern = "^cutadapt\\..*\\.json$", full.names = TRUE)

primer_stats <- list()
for (json_file in primer_jsons) {
  sample_name <- sub("^cutadapt\\.", "", sub("\\.json$", "", basename(json_file)))
  stats <- parse_cutadapt_json(json_file)
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
qc_df$adapter_r1_detected <- NA
qc_df$adapter_r1_pct <- NA
qc_df$adapter_r2_detected <- NA
qc_df$adapter_r2_pct <- NA
qc_df$adapter_passing <- NA
for (sname in sample_names) {
  if (sname %in% names(adapter_stats)) {
    qc_df[sname, "adapter_input"] <- adapter_stats[[sname]]$input
    qc_df[sname, "adapter_r1_detected"] <- adapter_stats[[sname]]$r1_with_adapter
    qc_df[sname, "adapter_r1_pct"] <- adapter_stats[[sname]]$pct_r1
    qc_df[sname, "adapter_r2_detected"] <- adapter_stats[[sname]]$r2_with_adapter
    qc_df[sname, "adapter_r2_pct"] <- adapter_stats[[sname]]$pct_r2
    qc_df[sname, "adapter_passing"] <- adapter_stats[[sname]]$output
  }
}

# Add primer stats
qc_df$primer_input <- NA
qc_df$primer_r1_detected <- NA
qc_df$primer_r1_pct <- NA
qc_df$primer_r2_detected <- NA
qc_df$primer_r2_pct <- NA
qc_df$primer_passing <- NA
for (sname in sample_names) {
  if (sname %in% names(primer_stats)) {
    qc_df[sname, "primer_input"] <- primer_stats[[sname]]$input
    qc_df[sname, "primer_r1_detected"] <- primer_stats[[sname]]$r1_with_adapter
    qc_df[sname, "primer_r1_pct"] <- primer_stats[[sname]]$pct_r1
    qc_df[sname, "primer_r2_detected"] <- primer_stats[[sname]]$r2_with_adapter
    qc_df[sname, "primer_r2_pct"] <- primer_stats[[sname]]$pct_r2
    qc_df[sname, "primer_passing"] <- primer_stats[[sname]]$output
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
cat("# COLUMN DESCRIPTIONS:\n")
cat("#   adapter_r*_pct:   % of R reads with Illumina adapter detected\n")
cat("#   primer_r*_pct:    % of R reads with amplicon primer detected\n")
cat("#   dada2_*:          DADA2 denoising results\n")
cat("#   merge_pct:        % of reads that successfully merged\n")
cat("#   chimera_pct:      % of merged reads retained after chimera removal\n")
cat("#\n")
cat("# EXPECTED VALUES (approximate):\n")
cat("#   adapter_r1/r2_pct:    typically 85-98% (low = adapter contamination or poor quality)\n")
cat("#   primer_r1/r2_pct:     typically 70-95% (low = wrong primer or orientation issues)\n")
cat("#                         For 18S-V4, expect R2 primer %% to be lower (50-70%)\n")
cat("#   dada2_filter_pct:     typically 70-85% (< 50% suggests quality problems)\n")
cat("#   merge_pct:            16S: 75-95%, ITS: 60-85%, 18S-V4: 30-70% (lower R2 quality)\n")
cat("#   chimera_pct:          typically 85-98% remaining (> 15% chimeric = concern)\n")
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
    p <- plotQualityProfile(primer_fqs_fwd[1:n_samples])
    ggsave(file.path(qc_outdir, "figures", "quality_profile_fwd.pdf"), plot = p, width = 12, height = 7)
    cat("  quality_profile_fwd.pdf\n")
  }, error = function(e) {
    cat("  Warning: could not generate quality_profile_fwd.pdf:", e$message, "\n")
  })

  tryCatch({
    p <- plotQualityProfile(primer_fqs_rev[1:n_samples])
    ggsave(file.path(qc_outdir, "figures", "quality_profile_rev.pdf"), plot = p, width = 12, height = 7)
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
  p <- plotErrors(errF, nominalQ = TRUE)
  ggsave(file.path(qc_outdir, "figures", "error_model_fwd.pdf"), plot = p, width = 12, height = 7)
  cat("  error_model_fwd.pdf\n")
}, error = function(e) {
  cat("  Warning: could not generate error_model_fwd.pdf:", e$message, "\n")
})

tryCatch({
  p <- plotErrors(errR, nominalQ = TRUE)
  ggsave(file.path(qc_outdir, "figures", "error_model_rev.pdf"), plot = p, width = 12, height = 7)
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
  library(ggplot2)

  # Calculate mean reads at each step
  steps <- c("adapter_input", "adapter_passing", "primer_passing", "dada2_input", "filtered", "merged", "nonchim")
  step_data <- data.frame()

  for (step in steps) {
    if (step %in% colnames(qc_df)) {
      means <- mean(qc_df[[step]], na.rm = TRUE)
    } else {
      means <- NA
    }
    step_data <- rbind(step_data, data.frame(step = step, mean_reads = means))
  }

  # Remove NAs
  step_data <- step_data[!is.na(step_data$mean_reads), ]

  # Create ggplot with better spacing
  p <- ggplot(step_data, aes(x = factor(step, levels = step), y = mean_reads, fill = step)) +
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
}, error = function(e) {
  cat("  Warning: could not generate read_tracking.pdf:", e$message, "\n")
})

cat("\nQC figures generated in:", file.path(qc_outdir, "figures"), "\n")
cat("\nDone.\n")
