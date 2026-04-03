library(dplyr)
# need to run "sort -h templistout.txt" numerically by first number
# Rscript output_asv.R seqtab_nochim.rds 
finalASVlista <- readLines("templistout.txt")
finalASVlist <- strsplit(finalASVlista, " ")
args = commandArgs(trailingOnly=TRUE)
ASVtab <- readRDS(args[1])
ASVtab = as.data.frame(t(ASVtab))

ASVtab2 = ASVtab
ASVtab2$asvgroup = 0

# starting at last item in list make all rows = group
endl <- length(finalASVlist)
# for each line in templistout
for (x in rev(seq_along(1:endl))) {
# get line
  these <- as.numeric(finalASVlist[[x]])
# for each item in list
  for (i in these){
  current <- ASVtab2[i, "asvgroup"]
  new = as.numeric(min(finalASVlist[[x]][1], current[current >0]))
  # this row switch to group (or stay same)
  ASVtab2[i,]$asvgroup <- new
  # any row in same group as this one should do same
  if(nrow(ASVtab2[ASVtab2$asvgroup == i,]>0)){
  ASVtab2[ASVtab2$asvgroup == i,]$asvgroup <- new
  }
  }
}

# split off rows that are not grouped
non_group <- subset(ASVtab2, asvgroup == 0)
yes_group <-  subset(ASVtab2, asvgroup != 0)

# grouped <- yes_group %>% group_by(asvgroup) %>% summarize_each(list(sum))
grouped <- yes_group %>% group_by(asvgroup) %>% summarize(across(where(is.numeric), sum), .groups = 'drop')

grouped <- as.data.frame(grouped)
for (x in 1:nrow(grouped)){
y = as.numeric(grouped[x,1])
row.names(grouped)[x] = row.names(ASVtab2[y,])
}

# recombine grouped and non-grouped asvs
grouped = grouped[,-c(1)]
end <- ncol(non_group)-1
non_group = non_group[,1:end]
both = rbind(grouped,non_group)

# save output
fname <- tools::file_path_sans_ext(args[1])
fname <- paste0(fname,"_merged_asvtable.rds")
saveRDS(both, file = fname)

# This removes asvs present in one sample only
# row_counts <- apply(both, 1, function(x) sum(x > 0))
# no_singles <- both[row_counts > 1,]
