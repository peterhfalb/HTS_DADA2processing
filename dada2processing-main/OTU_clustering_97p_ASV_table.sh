# install vsearch
module load conda
module load R/4.4.0-openblas-rocky8
source /common/software/install/migrated/anaconda/python3-2023.03-libmamba/etc/profile.d/conda.sh
conda create --name vsearch -c conda-forge -c bioconda bioconda::vsearch
conda activate vsearch

getSeqs.R
########################################
# gets sequences as a fasta from .rds file
# arg[1] = seqtab_nochim.R
args = commandArgs(trailingOnly=TRUE)
seqtab.nochim <- readRDS(args[1])
dada2::uniquesToFasta(seqtab.nochim, fout = "sequences.fasta")
########################################

otuSeqsID.R
########################################
# get OTU ID and sequence matching ID
# arg[1] = all.otutab.txt
args = commandArgs(trailingOnly=TRUE)
all.otutab <- read.delim(args[1])
library(tidyr)
library(dplyr)

writeLines(paste(c("X.OTU.ID","variable","value"), collapse = "\t"), "OTUID_SEQID.txt")
getLongform <- function(x){
data <- reshape2::melt(all.otutab[,c(1,x)], id.vars = c("X.OTU.ID"))
data[data==0] <- NA
data2<-data[complete.cases(data),]
write.table(data2,file="OTUID_SEQID.txt",append=TRUE, quote = FALSE, col.names = FALSE, row.names=FALSE, sep = "\t")
}
out <- lapply(2:ncol(all.otutab),getLongform)

########################################

OTUmerge.R
########################################
# arg1 = "seqtab_nochim.rds"
# arg2 = "OTUID_SEQID.txt"
args = commandArgs(trailingOnly=TRUE)
seqtab <- readRDS(args[1])
seqtab = t(seqtab)
OTUID <- read.delim(args[2])
seqtab.OTU <- data.frame(OTUID$X.OTU.ID,seqtab)
colnames(seqtab.OTU)[1] <- "OTU.ID"
out <- setDT(seqtab.OTU)[, lapply(.SD, sum), keyby = OTU.ID]
write.table(out, file = "OTU97percentID.txt", sep = "\t", quote = FALSE, row.names = FALSE)
########################################

taxaFromFasta.R
########################################
#' get taxonomy prediction of a fasta file
#' arg[1] is fasta file name
#' arg[2] is one of [16S|ITS|18S]
#' 
#' @export

args = commandArgs(trailingOnly=TRUE)

seqs <- Biostrings::readDNAStringSet(args[1])
#16S
ref16S <- "/home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa"
#18S
ref18S <- "/home/umii/public/dada2_taxonomy_references/maarjam_dada2.txt"
#ITS
ITSref <- "/home/umii/public/dada2_taxonomy_references/sh_general_release_dynamic_all_04.04.2024.fasta"

if (args[2]=="16S") {
taxa <- dada2::assignTaxonomy(seqs, ref16S, multithread = TRUE, outputBootstraps = TRUE)
} else if(args[2]=="18S") {
taxa <- dada2::assignTaxonomy(seqs, ref18S, multithread = TRUE, outputBootstraps = TRUE)
} else if(args[2]=="ITS"){
taxa <- dada2::assignTaxonomy(seqs, ITSref, multithread = TRUE, outputBootstraps = TRUE)
}
taxout <- taxa$tax
bootout <- taxa$boot
write.table(taxout, file = "taxaID.txt", sep = "\t", quote = FALSE)
write.table(bootout, file = "taxaIDbootstrap.txt", sep = "\t", quote = FALSE)
########################################
OTU97percentID <- read.delim("OTU97percentID.txt")
taxaID <- read.delim("taxaID.txt")
taxaID$OTU.ID <- paste0("OTU_",1:nrow(OTU97percentID))
both <- merge(OTU97percentID,taxaID, by = "OTU.ID")
write.table(both, file = "OTU97percentID_taxa.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# reordering by sum of row
OTU97percentID_taxa <- read.delim("OTU97percentID_taxa.txt")
OTU97percentID_taxa$SUMS <- rowSums(OTU97percentID_taxa[,c(2:224)])
sortedBoth <- OTU97percentID_taxa[order(OTU97percentID_taxa$SUMS, decreasing = TRUE),]
sortedBoth = sortedBoth[,-c(ncol(sortedBoth))]
write.table(sortedBoth, file = "OTU97percentID_taxa_sorted.txt", sep = "\t", quote = FALSE, row.names = FALSE)

########################################
OTU97percentID <- read.delim("OTU97percentID.txt")
blastN_out <- read.delim("blastN_out.tab", header=FALSE)
blastN_out <- blastN_out %>% tidyr::separate(V1, into = c("OTU.ID", "size"), sep = ";")
both <- merge(OTU97percentID,blastN_out, by = "OTU.ID", all = TRUE)
write.table(both, file = "OTU97percentID_blast.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# reordering by sum of row
both$SUMS <- rowSums(both[,c(2:224)])
sortedBoth <- both[order(both$SUMS, decreasing = TRUE),]
sortedBoth = sortedBoth[,-c(ncol(sortedBoth))]
write.table(sortedBoth, file = "OTU97percentID_blast_sorted.txt", sep = "\t", quote = FALSE, row.names = FALSE)

Running pipeline
########################################
module load conda
module load R/4.4.0-openblas-rocky8
source /common/software/install/migrated/anaconda/python3-2023.03-libmamba/etc/profile.d/conda.sh
conda activate vsearch

Rscript /home/umii/goul0109/scripts/getSeqs.R "seqtab_nochim.rds"
vsearch --cluster_size sequences.fasta --threads 8 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
Rscript /home/umii/goul0109/scripts/otuSeqsID.R "all.otutab.txt"
Rscript /home/umii/goul0109/scripts/OTUmerge.R "seqtab_nochim.rds" "OTUID_SEQID.txt"
# here you need to specify the [ITS/18S/16S] reference type
Rscript /home/umii/goul0109/scripts/taxaFromFasta.R "all.otus.fasta" "18S"

# important output files
# OTU97percentID.txt: OTU table with Samples x OTUs and cells are counts
# all.otus.fasta: contains all OTU centroid sequences (representative)
# taxaID.txt: taxonomic id of representative sequences 

#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=100g
#SBATCH --tmp=100g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

module load conda
module load R/4.4.0-openblas-rocky8
source /common/software/install/migrated/anaconda/python3-2023.03-libmamba/etc/profile.d/conda.sh
conda activate vsearch
# 1
cd /scratch.global/goul0109/kennedy_All_Dada2_Output/YucatanSoilsITS_DADA2_OUTPUT/
Rscript /home/umii/goul0109/scripts/getSeqs.R "seqtab_nochim.rds"
vsearch --cluster_size sequences.fasta --threads 24 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
Rscript /home/umii/goul0109/scripts/otuSeqsID.R "all.otutab.txt"
Rscript /home/umii/goul0109/scripts/OTUmerge.R "seqtab_nochim.rds" "OTUID_SEQID.txt"
Rscript /home/umii/goul0109/scripts/taxaFromFasta.R "all.otus.fasta" "ITS"
# 2
cd ../YucSoil16S_DADA2_OUTPUT/
Rscript /home/umii/goul0109/scripts/getSeqs.R "seqtab_nochim.rds"
vsearch --cluster_size sequences.fasta --threads 24 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
Rscript /home/umii/goul0109/scripts/otuSeqsID.R "all.otutab.txt"
Rscript /home/umii/goul0109/scripts/OTUmerge.R "seqtab_nochim.rds" "OTUID_SEQID.txt"
Rscript /home/umii/goul0109/scripts/taxaFromFasta.R "all.otus.fasta" "16S"
# 3
cd ../MycoSoil18S_DADA2_OUTPUT/
Rscript /home/umii/goul0109/scripts/getSeqs.R "seqtab_nochim.rds"
vsearch --cluster_size sequences.fasta --threads 24 --id 0.97 --strand plus --sizein --sizeout --fasta_width 0 --uc all.clustered.uc --relabel OTU_ --centroids all.otus.fasta --otutabout all.otutab.txt
Rscript /home/umii/goul0109/scripts/otuSeqsID.R "all.otutab.txt"
Rscript /home/umii/goul0109/scripts/OTUmerge.R "seqtab_nochim.rds" "OTUID_SEQID.txt"
Rscript /home/umii/goul0109/scripts/taxaFromFasta.R "all.otus.fasta" "18S"
# 4
