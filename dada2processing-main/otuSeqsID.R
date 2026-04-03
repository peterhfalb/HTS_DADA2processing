#' otuSeqsID
#' This script gets otu and sequence Id in same line from vsearch output
#' usage Rscript otuSeqsID.R "all.otutab.txt"
#'
#' @export

args = commandArgs(trailingOnly=TRUE)

all.otutab <- read.delim(args[1])
library(tidyr)
library(dplyr)

all.otutab[all.otutab == 0] <- NA
melted <- all.otutab %>% pivot_longer(cols = !X.OTU.ID, names_to = c("variable"), values_to = c("value"), values_drop_na = TRUE)
gsub("sq","",melted$variable)
melted <- melted[order(as.numeric(melted$variable)),]
write.table(melted, file = "OTUID_SEQID.txt", sep = "\t", quote = FALSE)