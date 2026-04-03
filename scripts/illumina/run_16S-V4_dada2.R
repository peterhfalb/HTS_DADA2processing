library(dada2)
args = commandArgs(trailingOnly=TRUE)

# Get output directory from args[2], default to ../dada2output/ for backwards compatibility
if (length(args) >= 2 && nchar(args[2]) > 0) {
  output_dir <- args[2]
} else {
  output_dir <- "../dada2output"
}

# Ensure output directory exists
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

path <- (".")
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Quality filtering based on args[1]
quality <- if (length(args) >= 1) args[1] else "good"

# 16S filtering parameters
if (quality == "good") {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), minLen = 100, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=128)
} else if (quality == "bad") {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(4,6), minLen = 100, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=8)
} else {
  stop("Quality must be 'good' or 'bad'")
}

head(out)

# Dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# Error models
errF <- learnErrors(derep_forward, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=TRUE, randomize=TRUE)

# DADA2 inference
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE, pool="pseudo")

# Merge pairs
merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=20)

# Make sequence table
seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Check sequence lengths
lname <- nchar(colnames(seqtab.nochim))
summary(lname)

saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

# Summary statistics
getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=sample.names,
                          dada2_input=out[,1],
                          filtered=out[,2],
                          dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN),
                          merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = file.path(output_dir, "sequence_process_summary.txt"), sep = "\t", quote=FALSE)

# Export sequences
seqtab.nochim <- readRDS(file.path(output_dir, "seqtab_nochim.rds"))
uniquesToFasta(seqtab.nochim, fout = file.path(output_dir, "sequences.fasta"))

# Taxonomy assignment using SILVA (Kennedy lab shared folder)
taxasilva <- assignTaxonomy(seqtab.nochim,
                            "/panfs/jay/groups/4/kennedyp/shared/taxonomy/silva_nr99_v138.1_train_set.fa.gz",
                            multithread=TRUE,
                            outputBootstraps = TRUE)
taxout <- taxasilva$tax
bootout <- taxasilva$boot
saveRDS(taxout, file = file.path(output_dir, "taxIDsilva.rds"))
saveRDS(bootout, file = file.path(output_dir, "taxIDsilva_bootstrap.rds"))

# Combine sequences with taxonomy
both1 <- cbind(t(seqtab.nochim), taxasilva$tax, taxasilva$boot)
both2 <- cbind(t(seqtab.nochim), taxasilva$tax)
write.table(both1, file = file.path(output_dir, "16S-V4_combined_sequences_taxa_silva_boot.txt"), sep = "\t", quote = FALSE, col.names=NA)
write.table(both2, file = file.path(output_dir, "16S-V4_combined_sequences_taxa_silva.txt"), sep = "\t", quote = FALSE, col.names=NA)

cat("16S-V4 DADA2 processing complete.\n")
cat("Output files written to:", output_dir, "\n")
