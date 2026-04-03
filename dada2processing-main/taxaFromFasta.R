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