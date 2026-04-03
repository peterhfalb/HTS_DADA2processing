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

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,160), maxN=0, maxEE=c(2,2), minLen = 100, truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# error models
errF <- learnErrors(derep_forward, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=TRUE, randomize=TRUE)

omegaset = paste0("1e-",args[1])
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE, pool="pseudo", OMEGA_A=omegaset)
dadaRs <- dada(derep_reverse, err=errR, multithread=TRUE, pool="pseudo", OMEGA_A=omegaset)

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=20)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "/seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

lname <- nchar(colnames(seqtab.nochim))
summary(lname)

filename1 = paste0("seqtab_nochim_omega",args[1],".rds")
saveRDS(seqtab.nochim, filename1)

  # set a little function
getN <- function(x) sum(getUniques(x))

  # making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN
               ),nonchim=rowSums(seqtab.nochim))
filename2 = paste0("sequence_process_summary_omega",args[1],".txt")
write.table(summary_tab, file = filename2, sep = "\t", quote=FALSE)

seqtab.nochim <- readRDS(filename1)

filename3 = paste0("sequences_omega",args[1],".fasta")
uniquesToFasta(seqtab.nochim, fout = filename3)

#TAXONOMY
taxasilva <- assignTaxonomy(seqtab.nochim, "/home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa", multithread=TRUE, outputBootstraps = TRUE)
taxout <- taxasilva$tax
bootout <- taxasilva$boot

filename4 = paste0("taxIDsilva_omega",args[1],".rds")
filename5 = paste0("taxIDsilva_bootstrap_omega",args[1],".rds")
saveRDS(taxout, file = filename4)
saveRDS(bootout, file = filename5)

both1 <- cbind(t(seqtab.nochim),taxasilva$tax, taxasilva$boot)
both2 <- cbind(t(seqtab.nochim),taxasilva$tax)

filename6 = paste0("16S_omega",args[1],"_combined_sequences_taxa_silva_boot.txt")
filename7 = paste0("16S_omega",args[1],"_combined_sequences_taxa_silva.txt")
write.table(both1, file = filename6, sep = "\t", quote = FALSE, col.names=NA)
write.table(both2, file = filename7, sep = "\t", quote = FALSE, col.names=NA)