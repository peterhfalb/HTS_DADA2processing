library(dada2)
args = commandArgs(trailingOnly=TRUE)

if (length(args) >= 2 && nchar(args[2]) > 0) {
  output_dir <- args[2]
} else {
  output_dir <- "../dada2output"
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

path <- (".")
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

quality <- if (length(args) >= 1) args[1] else "good"

# 18S-AMF filtering
if (quality == "good") {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=5, minLen = 100, maxEE=c(2,4), matchIDs=TRUE, rm.phix=TRUE, compress=TRUE, multithread=128)
} else if (quality == "bad") {
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=5, minLen = 100, maxEE=c(4,6), matchIDs=TRUE, rm.phix=TRUE, compress=TRUE, multithread=8)
} else {
  stop("Quality must be 'good' or 'bad'")
}

head(out)

derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

errF <- learnErrors(derep_forward, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=TRUE, randomize=TRUE)

dadaFs <- dada(derep_forward, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=10)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

getN <- function(x) sum(getUniques(x))
summary_tab <- data.frame(row.names=sample.names,
                          dada2_input=out[,1],
                          filtered=out[,2],
                          dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN),
                          merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = file.path(output_dir, "sequence_process_summary.txt"), sep = "\t", quote=FALSE)

seqtab.nochim <- readRDS(file.path(output_dir, "seqtab_nochim.rds"))
uniquesToFasta(seqtab.nochim, fout = file.path(output_dir, "sequences.fasta"))

# Taxonomy assignment using MaarjAM for AMF (Kennedy lab shared folder)
taxamaarjam <- assignTaxonomy(seqtab.nochim,
                              "/panfs/jay/groups/4/kennedyp/shared/taxonomy/MaarjAM_dada2_v.240.fa.gz",
                              multithread=TRUE,
                              outputBootstraps = TRUE)
taxout <- taxamaarjam$tax
bootout <- taxamaarjam$boot
saveRDS(taxout, file = file.path(output_dir, "taxIDmaarjam.rds"))
saveRDS(bootout, file = file.path(output_dir, "taxIDmaarjam_bootstrap.rds"))

both1 <- cbind(t(seqtab.nochim), taxamaarjam$tax, taxamaarjam$boot)
both2 <- cbind(t(seqtab.nochim), taxamaarjam$tax)
write.table(both1, file = file.path(output_dir, "18S-AMF_combined_sequences_taxa_maarjam_boot.txt"), sep = "\t", quote = FALSE, col.names=NA)
write.table(both2, file = file.path(output_dir, "18S-AMF_combined_sequences_taxa_maarjam.txt"), sep = "\t", quote = FALSE, col.names=NA)

cat("18S-AMF DADA2 processing complete.\n")
cat("Output files written to:", output_dir, "\n")
