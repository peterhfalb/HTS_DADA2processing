#!/usr/bin/env Rscript
# asv_amf_filter.R
# Applies the AMF SILVA screen to ASVs (18S-AMF only), before ASV taxonomy assignment.
# The keep/discard decisions are made by amf_silva_filter.R, as for OTUs; this script
# converts the DADA2 sequence table to and from the inputs and outputs of that step.
#
# Usage: Rscript asv_amf_filter.R prepare <seqtab_rds> <work_dir>
#        Rscript asv_amf_filter.R finalize <seqtab_rds> <work_dir> <output_rds> <project>

library(readr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || !args[1] %in% c("prepare", "finalize")) {
  stop("Usage: Rscript asv_amf_filter.R prepare <seqtab_rds> <work_dir>\n",
       "       Rscript asv_amf_filter.R finalize <seqtab_rds> <work_dir> <output_rds> <project>")
}

MODE     <- args[1]
WORK_DIR <- args[3]
seqtab   <- readRDS(args[2])

asv_seqs <- colnames(seqtab)
asv_ids  <- sprintf("ASV%06d", seq_along(asv_seqs))

fasta_path  <- file.path(WORK_DIR, "asv_sequences.fasta")
counts_path <- file.path(WORK_DIR, "asv_counts_with_ids.txt")
map_path    <- file.path(WORK_DIR, "asv_empty_map.csv")

if (MODE == "prepare") {
  dir.create(WORK_DIR, showWarnings = FALSE, recursive = TRUE)
  writeLines(as.vector(rbind(paste0(">", asv_ids), asv_seqs)), fasta_path)
  counts <- data.frame(ASV_ID = asv_ids, t(seqtab), check.names = FALSE)
  write.table(counts, counts_path, sep = "\t", quote = FALSE, row.names = FALSE)
  # Sample columns here already hold the original sample names, so nothing to map
  writeLines("BioInf_Sample,Sample_name", map_path)
  cat("Prepared", length(asv_ids), "ASVs for the AMF SILVA screen\n")
  quit(status = 0)
}

if (length(args) < 5) {
  stop("finalize needs <output_rds> and <project>")
}
OUTPUT_RDS <- args[4]
PROJECT    <- args[5]

# Label the decision files with ASV IDs and sequences, since the ASV tables are keyed by sequence
add_sequences <- function(path) {
  d <- read_tsv(path, show_col_types = FALSE)
  names(d)[names(d) == "OTU_ID"] <- "ASV_ID"
  d$ASV_sequence <- asv_seqs[match(d$ASV_ID, asv_ids)]
  write_tsv(d, path, na = "")
  d
}
decisions <- add_sequences(file.path(WORK_DIR, "amf_filter_otu_decisions.tsv"))
add_sequences(file.path(WORK_DIR, "amf_filter_discarded.tsv"))
file.rename(file.path(WORK_DIR, "amf_filter_otu_decisions.tsv"),
            file.path(WORK_DIR, "amf_filter_asv_decisions.tsv"))
file.rename(file.path(WORK_DIR, "Centroid_AMF_filtered.fas"),
            file.path(WORK_DIR, "ASVs_AMF_filtered.fasta"))

kept <- decisions$ASV_sequence[decisions$decision == "keep"]
saveRDS(seqtab[, asv_seqs %in% kept, drop = FALSE], OUTPUT_RDS)
cat("Filtered sequence table:", sum(asv_seqs %in% kept), "of", length(asv_seqs), "ASVs\n")

# Intermediate files, all reproducible from the outputs above
unlink(c(counts_path, map_path, fasta_path, file.path(WORK_DIR, paste0(PROJECT, "_AMF_filtered.txt"))))
