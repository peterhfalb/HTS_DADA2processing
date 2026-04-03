#' getSeqs
#' This script gets sequences as a fasta from .rds file
#' usage Rscript getSeqs.R "seqtab_nochim.rds"
#'
#' @export

args = commandArgs(trailingOnly=TRUE)
library(dada2)
seqtab.nochim <- readRDS(args[1])
uniquesToFasta(seqtab.nochim, fout = "sequences.fasta")