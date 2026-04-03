library(dada2)
args = commandArgs(trailingOnly=TRUE)
path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.paired.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.paired.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

quality=args[1]
# good quality
if (quality == "good"){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=5, minLen = 100, maxEE=c(2,4), matchIDs=TRUE, maxN = 0, rm.phix=TRUE, multithread=TRUE, verbose = TRUE)
}
# bad quality
if (quality == "bad"){
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncQ=5, minLen = 100, maxEE=c(4,6), matchIDs=TRUE, maxN = 0, rm.phix=TRUE, multithread=TRUE, verbose = TRUE)
}

head(out)

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# get string of qualities
files <- list.files(pattern = "\\.fastq.gz$")
fastq <- readFastq(files[1])
quals <- quality(fastq)
# Convert to a character matrix of ASCII-encoded quality characters
qual_chars <- as(quals, "matrix")
binned_Qs <- unique(as.vector((qual_chars)))
BQEF <- makeBinnedQualErrfun(binned_Qs)

# error models
errF <- learnErrors(derep_forward, errorEstimationFunction=BQEF, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(derep_reverse, errorEstimationFunction=BQEF, multithread=TRUE, randomize=TRUE)

dadaFs <- dada(derep_forward, err=errF, multithread=TRUE, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, minOverlap=10)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "../dada2output/seqtab.rds")

# Assumes seqtab is your sequence table of merged sequences
#MINLEN <- 400
#MAXLEN <- 600
#seqlens <- nchar(getSequences(seqtab))
#seqtab.filt <- seqtab[,seqlens >= MINLEN & seqlens <= MAXLEN]
#seqtab.nochim <- removeBimeraDenovo(seqtab.filt, method="consensus", multithread=TRUE, verbose=TRUE)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

lname <- nchar(colnames(seqtab.nochim))
summary(lname)

dim(seqtab.nochim)
saveRDS(seqtab.nochim, "../dada2output/seqtab_nochim.rds")
uniquesToFasta(seqtab.nochim, fout = "../dada2output/sequences.fasta")

  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "../dada2output/sequence_process_summary.txt", sep = "\t", quote=FALSE)

taxrefa <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
taxa <- assignTaxonomy(seqtab.nochim, taxrefa, tryRC = TRUE, taxLevels = c("Class", "Order", "Family", "Genus", "Species"), multithread = TRUE, outputBootstraps = TRUE)
taxout <- taxa$tax
bootout <- taxa$boot
saveRDS(taxout, file = "../dada2output/taxIDmaar.rds")
saveRDS(bootout, file = "../dada2output/taxIDmaar_bootstrap.rds")
both1 <- cbind(t(seqtab.nochim),taxout,bootout)
write.table(both1, file = "../dada2output/18S_combined_sequences_taxamaar_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxout)
write.table(both2, file = "../dada2output/18S_combined_sequences_taxamaar.txt", sep = "\t", quote = FALSE, col.names=NA)
