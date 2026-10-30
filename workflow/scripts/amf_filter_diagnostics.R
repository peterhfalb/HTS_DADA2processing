#!/usr/bin/env Rscript
# amf_filter_diagnostics.R
# Diagnostics for the 18S-AMF SILVA screen (called by bin/diagnose_amf_filter.sh).
#
# Usage:
#   Rscript amf_filter_diagnostics.R select <decisions_tsv> <curated_fasta> <diag_dir> <n_top>
#     Writes reason_summary.tsv and top_<reason>.fasta (most abundant OTUs per discard reason)
#   Rscript amf_filter_diagnostics.R report <decisions_tsv> <diag_dir> <target_regex>
#     Reads blast_<reason>_<db>.tsv (unfiltered BLAST of the top OTUs) and writes
#     top_<reason>_besthits.tsv with the best hit per OTU in each database
#
# BLAST columns: qseqid sseqid pident length qlen qcovhsp bitscore evalue stitle

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
MODE <- if (length(args) >= 1) args[1] else ""

DIAG_REASONS <- c("no_confident_hit", "non_AMF")

read_fasta_lengths <- function(path) {
  lines  <- readLines(path)
  hdr    <- grep("^>", lines)
  ends   <- c(hdr[-1] - 1, length(lines))
  ids    <- sub(";.*", "", sub("^>", "", lines[hdr]))
  seqs   <- vapply(seq_along(hdr), function(i) {
    if (ends[i] > hdr[i]) paste(lines[(hdr[i] + 1):ends[i]], collapse = "") else ""
  }, character(1))
  tibble(OTU_ID = ids, seq = seqs, length = nchar(seqs))
}

# ------------------------------------------------------------------------------
# select: per-reason summary and query fastas
# ------------------------------------------------------------------------------
if (MODE == "select") {
  if (length(args) < 5) stop("Usage: amf_filter_diagnostics.R select <decisions_tsv> <curated_fasta> <diag_dir> <n_top>")
  DECISIONS <- args[2]
  FASTA     <- args[3]
  DIAG_DIR  <- args[4]
  N_TOP     <- as.integer(args[5])

  decisions <- read_tsv(DECISIONS, show_col_types = FALSE, col_types = cols(OTU_ID = "c"))
  seqs      <- read_fasta_lengths(FASTA)
  decisions <- decisions %>% left_join(seqs, by = "OTU_ID")

  reads_all <- sum(decisions$total_reads)
  reason_summary <- decisions %>%
    group_by(reason) %>%
    summarise(
      n_otus      = n(),
      reads       = sum(total_reads),
      pct_reads   = round(100 * sum(total_reads) / reads_all, 1),
      len_min     = min(length, na.rm = TRUE),
      len_median  = median(length, na.rm = TRUE),
      len_max     = max(length, na.rm = TRUE),
      # identity/coverage of best hit, only defined for OTUs with a hit passing the cutoffs
      best_pident_median = suppressWarnings(median(best_pident, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    arrange(desc(reads))

  write_tsv(reason_summary, file.path(DIAG_DIR, "reason_summary.tsv"), na = "")
  cat("\n== OTUs and reads by AMF filter decision ==\n")
  print(as.data.frame(reason_summary), row.names = FALSE)

  # Length distribution per reason (coarse bins), to spot odd-length artefacts
  len_bins <- decisions %>%
    mutate(len_bin = cut(length, breaks = c(0, 200, 300, 400, 450, 500, 550, 600, Inf), right = FALSE)) %>%
    count(reason, len_bin) %>%
    arrange(len_bin) %>%
    pivot_wider(names_from = reason, values_from = n, values_fill = 0)
  write_tsv(len_bins, file.path(DIAG_DIR, "length_bins_by_reason.tsv"))
  cat("\n== OTU length distribution by decision ==\n")
  print(as.data.frame(len_bins), row.names = FALSE)

  for (r in DIAG_REASONS) {
    top <- decisions %>% filter(reason == r, !is.na(seq)) %>% arrange(desc(total_reads)) %>% head(N_TOP)
    out <- file.path(DIAG_DIR, paste0("top_", r, ".fasta"))
    if (nrow(top) > 0) {
      writeLines(as.vector(rbind(paste0(">", top$OTU_ID), top$seq)), out)
    } else if (file.exists(out)) {
      file.remove(out)
    }
    cat("\nTop", r, "OTUs selected for unfiltered BLAST:", nrow(top), "\n")
  }

# ------------------------------------------------------------------------------
# report: best hit per OTU per database
# ------------------------------------------------------------------------------
} else if (MODE == "report") {
  if (length(args) < 4) stop("Usage: amf_filter_diagnostics.R report <decisions_tsv> <diag_dir> <target_regex>")
  DECISIONS    <- args[2]
  DIAG_DIR     <- args[3]
  TARGET_REGEX <- args[4]

  decisions <- read_tsv(DECISIONS, show_col_types = FALSE, col_types = cols(OTU_ID = "c"))
  blast_cols <- c("qseqid", "sseqid", "pident", "length", "qlen", "qcovhsp", "bitscore", "evalue", "stitle")

  for (r in DIAG_REASONS) {
    files <- list.files(DIAG_DIR, pattern = paste0("^blast_", r, "_.*\\.tsv$"), full.names = TRUE)
    query <- file.path(DIAG_DIR, paste0("top_", r, ".fasta"))
    if (!file.exists(query)) next
    ids <- sub(";.*", "", sub("^>", "", grep("^>", readLines(query), value = TRUE)))

    best <- bind_rows(lapply(files, function(f) {
      db <- sub(paste0("^blast_", r, "_(.*)\\.tsv$"), "\\1", basename(f))
      if (file.info(f)$size == 0) return(NULL)
      read_tsv(f, col_names = blast_cols, quote = "", show_col_types = FALSE,
               col_types = cols(.default = "c", pident = "d", length = "d", qlen = "d",
                                qcovhsp = "d", bitscore = "d", evalue = "d")) %>%
        mutate(db = db)
    }))

    if (nrow(best) > 0) {
      # Same taxonomy parsing as amf_silva_filter.R (strip a leading accession if present)
      best <- best %>%
        mutate(OTU_ID   = sub(";.*", "", qseqid),
               taxonomy = sub("^[^;[:space:]]+[[:space:]]+(?=[^;]*;)", "", stitle, perl = TRUE),
               amf      = grepl(TARGET_REGEX, taxonomy, ignore.case = TRUE)) %>%
        group_by(OTU_ID, db) %>%
        arrange(desc(bitscore), desc(pident), .by_group = TRUE) %>%
        summarise(pident = first(pident), qcov = first(qcovhsp), amf = first(amf),
                  taxonomy = first(taxonomy), .groups = "drop")
    } else {
      best <- tibble(OTU_ID = character(), db = character(), pident = numeric(),
                     qcov = numeric(), amf = logical(), taxonomy = character())
    }
    dbs <- sort(unique(best$db))

    wide <- tibble(OTU_ID = ids) %>%
      left_join(select(decisions, OTU_ID, total_reads), by = "OTU_ID")
    for (d in dbs) {
      b <- best %>% filter(db == d) %>% select(OTU_ID, pident, qcov, amf, taxonomy)
      names(b)[-1] <- paste0(d, "_", names(b)[-1])
      wide <- wide %>% left_join(b, by = "OTU_ID")
    }
    write_tsv(wide, file.path(DIAG_DIR, paste0("top_", r, "_besthits.tsv")), na = "")

    cat("\n================================================================\n")
    cat("Top", length(ids), r, "OTUs: best unfiltered hit per database\n")
    cat("(pident/qcov, * = best hit matches '", TARGET_REGEX, "')\n", sep = "")
    cat("================================================================\n")
    compact <- wide %>% select(OTU_ID, total_reads)
    for (d in dbs) {
      p <- wide[[paste0(d, "_pident")]]; q <- wide[[paste0(d, "_qcov")]]; a <- wide[[paste0(d, "_amf")]]
      compact[[d]] <- ifelse(is.na(p), "-", sprintf("%.1f/%.0f%s", p, q, ifelse(a, "*", "")))
    }
    if ("SILVA" %in% dbs) {
      # Last two ranks of the best SILVA hit, for a quick look at what it is
      compact$SILVA_best <- vapply(strsplit(ifelse(is.na(wide$SILVA_taxonomy), "", wide$SILVA_taxonomy), ";"),
                                   function(x) paste(tail(x, 2), collapse = ";"), character(1))
    }
    print(as.data.frame(compact), row.names = FALSE)

    cat("\nPer database (of", length(ids), "OTUs):\n")
    for (d in dbs) {
      b <- best %>% filter(db == d)
      cat(sprintf("  %-10s hit: %3d | best hit AMF: %3d | pident >=97: %3d | 90-97: %3d | <90: %3d | qcov <80: %3d\n",
                  d, nrow(b), sum(b$amf), sum(b$pident >= 97), sum(b$pident >= 90 & b$pident < 97),
                  sum(b$pident < 90), sum(b$qcov < 80)))
    }
  }

  cat("\nHow to read this:\n")
  cat("  (a) Divergent AMF:  SILVA best hit AMF (*) at <90%, MaarjAM/Eukaryome >=95%\n")
  cat("      -> SILVA lacks close AMF references; add MaarjAM / use Eukaryome, or lower amf_min_id\n")
  cat("  (b) Coverage:       high pident but qcov <80, or odd lengths in length_bins_by_reason.tsv\n")
  cat("      -> chimeras/untrimmed ends; check upstream trimming or amf_min_qcov\n")
  cat("  (c) Non-18S:        no hit, or low pident everywhere\n")
  cat("      -> genuinely off-target; discarding them is correct\n")
  cat("  non_AMF OTUs with a high-identity AMF hit in MaarjAM/Eukaryome may be AMF losing the\n")
  cat("  majority vote to non-AMF SILVA references (check amf_top_pct / amf_min_frac)\n")

} else {
  stop("Usage: amf_filter_diagnostics.R select|report ...")
}
