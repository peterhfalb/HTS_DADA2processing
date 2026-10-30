#!/usr/bin/env Rscript
# amf_silva_filter.R
# Keeps only OTUs whose near-top SILVA BLAST hits are mostly AMF (18S-AMF only).
#
# Usage: Rscript amf_silva_filter.R <blast_tsv> <centroids_fasta> <otu_table> <map_file> <out_dir> \
#                                   <project> <top_pct> <min_frac> <target_regex>
#
# blast_tsv columns (blastn -outfmt "6 qseqid sseqid pident length qcovhsp bitscore evalue stitle"),
# already restricted to hits passing the identity and query coverage cutoffs.

library(dplyr)
library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
  stop("Usage: Rscript amf_silva_filter.R <blast_tsv> <centroids_fasta> <otu_table> <map_file> <out_dir> <project> <top_pct> <min_frac> <target_regex>")
}

BLAST_TSV    <- args[1]
FASTA        <- args[2]
OTU_TABLE    <- args[3]
MAP_FILE     <- args[4]
OUT_DIR      <- args[5]
PROJECT      <- args[6]
TOP_PCT      <- as.numeric(args[7])
MIN_FRAC     <- as.numeric(args[8])
TARGET_REGEX <- args[9]

cat("BLAST results:  ", BLAST_TSV, "\n")
cat("Centroids:      ", FASTA, "\n")
cat("OTU table:      ", OTU_TABLE, "\n")
cat("Near-top window:", TOP_PCT, "% of best bitscore\n")
cat("Min AMF frac:   ", MIN_FRAC, "\n")
cat("Target regex:   ", TARGET_REGEX, "\n\n")

# ------------------------------------------------------------------------------
# Load inputs
# ------------------------------------------------------------------------------

blast_cols <- c("qseqid", "sseqid", "pident", "length", "qcovhsp", "bitscore", "evalue", "stitle")
if (file.exists(BLAST_TSV) && file.info(BLAST_TSV)$size > 0) {
  blast <- read_tsv(BLAST_TSV, col_names = blast_cols, quote = "",
                    col_types = cols(.default = "c", pident = "d", length = "d",
                                     qcovhsp = "d", bitscore = "d", evalue = "d"))
} else {
  blast <- tibble(qseqid = character(), sseqid = character(), pident = numeric(),
                  length = numeric(), qcovhsp = numeric(), bitscore = numeric(),
                  evalue = numeric(), stitle = character())
}
blast <- blast %>% mutate(OTU_ID = sub(";.*", "", qseqid))
cat("BLAST hits loaded:", nrow(blast), "\n")

fasta_lines <- readLines(FASTA)
header_idx  <- grep("^>", fasta_lines)
fasta_ids   <- sub(";.*", "", sub("^>", "", fasta_lines[header_idx]))

otu <- read_tsv(OTU_TABLE, show_col_types = FALSE)
names(otu)[1] <- "OTU_ID"
sample_cols <- names(otu)[-1]
otu_reads <- tibble(OTU_ID = otu$OTU_ID, total_reads = rowSums(otu[, sample_cols, drop = FALSE]))

map <- read_csv(MAP_FILE, show_col_types = FALSE)

cat("OTUs in curated table:", nrow(otu), "\n\n")

# ------------------------------------------------------------------------------
# Decide per OTU
# ------------------------------------------------------------------------------

# SILVA titles are "<accession> <taxonomy path>" (species names contain spaces).
# Strip a leading accession token only if present (it has no ';'), keep the taxonomy part
blast <- blast %>%
  mutate(taxonomy = sub("^[^;[:space:]]+[[:space:]]+(?=[^;]*;)", "", stitle, perl = TRUE),
         is_target = grepl(TARGET_REGEX, taxonomy, ignore.case = TRUE))

hit_summary <- blast %>%
  group_by(OTU_ID) %>%
  mutate(best_bitscore = max(bitscore)) %>%
  filter(bitscore >= best_bitscore * (1 - TOP_PCT / 100)) %>%
  arrange(desc(bitscore), desc(pident), .by_group = TRUE) %>%
  summarise(
    best_hit      = first(sseqid),
    best_pident   = first(pident),
    best_qcov     = first(qcovhsp),
    best_bitscore = first(bitscore),
    best_taxonomy = first(taxonomy),
    n_top_hits    = n(),
    amf_frac      = mean(is_target),
    .groups = "drop"
  )

decisions <- otu_reads %>%
  left_join(hit_summary, by = "OTU_ID") %>%
  mutate(
    decision = ifelse(!is.na(amf_frac) & amf_frac >= MIN_FRAC, "keep", "discard"),
    reason   = case_when(
      is.na(amf_frac)       ~ "no_confident_hit",
      amf_frac >= MIN_FRAC  ~ "AMF",
      TRUE                  ~ "non_AMF"
    )
  ) %>%
  select(OTU_ID, total_reads, decision, reason, amf_frac, n_top_hits,
         best_hit, best_pident, best_qcov, best_bitscore, best_taxonomy)

kept_ids <- decisions$OTU_ID[decisions$decision == "keep"]

n_total     <- nrow(decisions)
n_kept      <- length(kept_ids)
n_non_amf   <- sum(decisions$reason == "non_AMF")
n_no_hit    <- sum(decisions$reason == "no_confident_hit")
reads_total <- sum(decisions$total_reads)
reads_kept  <- sum(decisions$total_reads[decisions$decision == "keep"])

cat("OTUs kept:               ", n_kept, "/", n_total, "\n")
cat("OTUs discarded (non-AMF):", n_non_amf, "\n")
cat("OTUs discarded (no hit): ", n_no_hit, "\n")
cat("Reads kept:              ", reads_kept, "/", reads_total, "\n\n")

if (n_kept == 0) {
  stop("No OTUs passed the AMF SILVA filter. Check the BLAST database, target regex (",
       TARGET_REGEX, ") and amf_filter_otu_decisions.tsv")
}

# ------------------------------------------------------------------------------
# Write outputs
# ------------------------------------------------------------------------------

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

write_tsv(decisions, file.path(OUT_DIR, "amf_filter_otu_decisions.tsv"), na = "")
write_tsv(filter(decisions, decision == "discard"), file.path(OUT_DIR, "amf_filter_discarded.tsv"), na = "")

# Filtered OTU table (same format as the mumu curated table)
otu_kept <- otu %>% filter(OTU_ID %in% kept_ids)
names(otu_kept)[1] <- names(read_tsv(OTU_TABLE, n_max = 0, show_col_types = FALSE))[1]
write.table(otu_kept, file.path(OUT_DIR, paste0(PROJECT, "_AMF_filtered.txt")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Filtered fasta, original header lines (with ;size=) preserved
seq_end  <- c(header_idx[-1] - 1, length(fasta_lines))
keep_rec <- which(fasta_ids %in% kept_ids)
out_lines <- unlist(lapply(keep_rec, function(i) fasta_lines[header_idx[i]:seq_end[i]]))
writeLines(out_lines, file.path(OUT_DIR, "Centroid_AMF_filtered.fas"))

# Per-sample read retention (sample columns may be S001-style IDs or original names)
per_sample <- tibble(
  sample      = sample_cols,
  reads_total = colSums(otu[, sample_cols, drop = FALSE]),
  reads_kept  = colSums(otu_kept[, sample_cols, drop = FALSE])
) %>%
  mutate(
    sample   = ifelse(sample %in% map$BioInf_Sample,
                      map$Sample_name[match(sample, map$BioInf_Sample)], sample),
    pct_kept = ifelse(reads_total > 0, round(100 * reads_kept / reads_total, 1), NA)
  )
write_tsv(per_sample, file.path(OUT_DIR, "amf_filter_per_sample.tsv"), na = "")

# Counts for the OTU QC summary (key=value, like the VSEARCH counts file)
pct <- function(x, n) if (n > 0) sprintf("%.1f", 100 * x / n) else "NA"
writeLines(c(
  paste0("AMF_INPUT_OTUS=", n_total),
  paste0("AMF_KEPT_OTUS=", n_kept),
  paste0("AMF_DISCARDED_NON_AMF=", n_non_amf),
  paste0("AMF_DISCARDED_NO_HIT=", n_no_hit),
  paste0("AMF_PCT_OTUS_DISCARDED=", pct(n_total - n_kept, n_total)),
  paste0("AMF_READS_TOTAL=", reads_total),
  paste0("AMF_READS_KEPT=", reads_kept),
  paste0("AMF_PCT_READS_RETAINED=", pct(reads_kept, reads_total))
), file.path(OUT_DIR, "amf_filter_counts.txt"))

cat("AMF filter outputs written to:", OUT_DIR, "\n")
