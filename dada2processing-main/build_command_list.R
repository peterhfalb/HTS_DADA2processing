#Rscript build_command_list.R seqtab_nochim.rds [integer errors allowed]
args = commandArgs(trailingOnly=TRUE)
ASVtab <- readRDS(args[1])
ASVtab = as.data.frame(t(ASVtab))
endnum = nrow(ASVtab)-1
# set default number of errors allowed
# wrong = 2
wrong = as.numeric(args[2])
for (i in (1:endnum)){
comlist <- paste0("Rscript /users/4/goul0109/dada2processing/clust_asv.R ",args[1]," ",i," ",wrong)
lapply(comlist, write, "runcommands.cmd", append=TRUE)
}