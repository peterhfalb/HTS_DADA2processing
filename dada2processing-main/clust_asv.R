# Rscript clust_asv.R seqtab_nochim.rds 
args = commandArgs(trailingOnly=TRUE)
#
ASVtab <- readRDS(args[1])
ASVtab = as.data.frame(t(ASVtab))

# set default number of errors allowed
wrong = as.numeric(args[3])

mat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
# function for the parallel loop
getalignrow <- function(x){
  seqs <- row.names(ASVtab)
  ref_seq = seqs[x]
  start = x+1
  query <- row.names(ASVtab[start:nrow(ASVtab),,drop=FALSE])
  #out <- pairwiseAlignment(query, ref_seq, substitutionMatrix = mat, type = "global", gapOpening = 5, gapExtension = 2, scoreOnly=TRUE)
  out <- pwalign::pairwiseAlignment(query, ref_seq, substitutionMatrix = mat, type = "global", gapOpening = 5, gapExtension = 2, scoreOnly=TRUE)
  subASV = ASVtab
  subASV$score = 0
  subASV[start:nrow(ASVtab),]$score = out

  max_score = nchar(row.names(ASVtab[1,,drop=FALSE]))
  min_score = nchar(row.names(ASVtab[1,,drop=FALSE])) - (4 * wrong)
  
  subASV$rownumber =1:nrow(subASV)
  
  subASV = subset(subASV, subASV$score >= min_score)
  c <- row.names(subASV[1,])
  d=wrong+1
  d <- stringdist::stringdist(ref_seq, c, method = "lv")
  print(x)
  if (d[1]<=wrong){
    mergerows <- c(x,subASV$rownumber)
    mergerows <- t(mergerows)
   return(mergerows)
  }
}
i = as.numeric(args[2])
alignlist <- getalignrow(i)
filename1 = paste0(args[1],"_templistout.txt")
write.table(alignlist, filename1, append=TRUE, row.names = FALSE, col.names = FALSE)