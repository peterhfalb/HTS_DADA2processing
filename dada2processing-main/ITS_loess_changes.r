library(dada2)
library(dplyr)

path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# ITS
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)

# dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

# error models
errF <- learnErrors(derep_forward, multithread=24, randomize=TRUE)
errR <- learnErrors(derep_reverse, multithread=24, randomize=TRUE)

ferrplot <- plotErrors(errF, nominalQ=TRUE)
rerrplot <- plotErrors(errR, nominalQ=TRUE)
ggplot2::ggsave(ferrplot, file = "Forward_error_plot.png", dpi = 600, height = 8, width = 8, units = "in")
ggplot2::ggsave(rerrplot, file = "Reverse_error_plot.png", dpi = 600, height = 8, width = 8, units = "in")
################################

# Source: https://github.com/benjjneb/dada2/issues/1307#issuecomment-957680971
loessErrfun_mod4 <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # jonalim's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),degree = 1, span = 0.95)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF_4 <- learnErrors(
  filtFs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

errR_4 <- learnErrors(
  filtRs,
  multithread = TRUE,
  nbases = 1e8,
  errorEstimationFunction = loessErrfun_mod4,
  verbose = TRUE
)

##########################

errFplot <- plotErrors(errF_4, nominalQ=TRUE)
errRplot <- plotErrors(errR_4, nominalQ=TRUE)
ggplot2::ggsave(errFplot, file = "Forward_error_plot_loess.png", dpi = 600, height = 8, width = 8, units = "in")
ggplot2::ggsave(errRplot, file = "Reverse_error_plot_loess.png", dpi = 600, height = 8, width = 8, units = "in")

dadaFs <- dada(derep_forward, err=errF_4, multithread=8, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR_4, multithread=8, pool="pseudo")

fderrplot <- plotErrors(dadaFs, nominalQ=TRUE)
rderrplot <- plotErrors(dadaRs, nominalQ=TRUE)

ggplot2::ggsave(fderrplot, file = "Forward_post_error_plot_loess.png", dpi = 600, height = 8, width = 8, units = "in")
ggplot2::ggsave(rderrplot, file = "Reverse_post_error_plot_loess.png", dpi = 600, height = 8, width = 8, units = "in")

# merging pairs
merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=50)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "../dada2output/seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "../dada2output/seqtab_nochim.rds")

# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
               filtered=out[,2], dada_f=sapply(dadaFs, getN),
               dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),nonchim=rowSums(seqtab.nochim))
write.table(summary_tab, file = "sequence_process_summary.txt", sep = "\t", quote=FALSE)

#####
ref <- "/panfs/jay/groups/4/kennedyp/shared/taxonomy/sh_general_release_dynamic_s_all_25.07.2023_dev.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ref, multithread = 8, outputBootstraps = TRUE)
saveRDS(taxa$tax, file = "../dada2output/taxa.rds")
saveRDS(taxa$boot, file = "../dada2output/taxa_bootstrap.rds")
#
both1 <- cbind(t(seqtab.nochim),taxa$tax, taxa$boot)
write.table(both1, file = "../dada2output/ITS_combined_sequences_taxa_bootstrap.txt", sep = "\t", quote = FALSE, col.names=NA)
both2 <- cbind(t(seqtab.nochim),taxa$tax)
write.table(both2, file = "../dada2output/ITS_combined_sequences_taxa.txt", sep = "\t", quote = FALSE, col.names=NA)
