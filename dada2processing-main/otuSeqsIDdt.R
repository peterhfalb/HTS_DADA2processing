library(data.table)

all.otutab <- fread(args[1])  # Faster and uses less memory
fwrite(data.table(X.OTU.ID = "X.OTU.ID", variable = "variable", value = "value"), 
       file = "OTUID_SEQID.txt", sep = "\t", col.names = FALSE)

for (i in 2:ncol(all.otutab)) {
  dt <- all.otutab[, .(X.OTU.ID, val = .SD[[1]]), .SDcols = i]
  setnames(dt, "val", names(all.otutab)[i])
  dt[dt[[2]] == 0, (2) := NA]
  dt <- na.omit(melt(dt, id.vars = "X.OTU.ID", variable.name = "variable"))
  fwrite(dt, file = "OTUID_SEQID.txt", sep = "\t", append = TRUE, col.names = FALSE)
  rm(dt); gc()
}