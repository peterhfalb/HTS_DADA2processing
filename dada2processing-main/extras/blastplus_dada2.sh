#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=128
#SBATCH --mem=499g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

# https://archive.jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/taxonomy-guide/

cd /working/Directory/
module load ncbi_blast+
blastn -db core_nt -query sequences.fasta -outfmt '6 qseqid sseqid evalue bitscore sgi sacc staxids sscinames scomnames stitle sblastnames sskingdoms'  -max_target_seqs 1  -max_hsps 1 -num_threads 128 > blastN_out.tab
cat /home/umii/goul0109/scripts/standard_blast_fields.tsv blastN_out.tab > t && mv t blastN_out_tabbed.tab

#!/bin/bash -l
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=200g
#SBATCH --tmp=200g
#SBATCH --mail-type=ALL
#SBATCH --mail-user=goul0109@umn.edu

module load R/4.4.0-openblas-rocky8

# get accession numbers
awk '{print $7}' blastN_out.tab > asc.txt
# remove rest of line after ; in any row with a ;
sed -i 's/;.*//g' asc.txt

# get taxonomies from accession numbers
/home/umii/goul0109/bbmap/taxonomy.sh tree=/home/umii/goul0109/tree.taxtree.gz in=asc.txt out=all_taxa.txt

# remove first 4 lines
sed -i '1,4d' all_taxa.txt

# reformat output into a table
mkdir temp_split
awk -v RS= '{print > ("temp_split/" NR ".txt")}' all_taxa.txt
cd temp_split
for file in *.txt; do 
base=$(basename "$file" '.txt')
cat $file | 
tac | 
grep -w 'kingdom\|phylum\|class\|order\|family\|genus\|species' | 
sed 's/\t/_/2g;P;D;' | 
sed "s/$/\t $base/g" >> ../long.txt; done
cd ..
rm -r temp_split

Rscript /home/umii/goul0109/scripts/long2wide.R long.txt

# long2wide.R
#########################################################################
#' long2wide
#' This script takes long format ncbiout to wide format table
#' usage Rscript long2wide.R "long.txt"
#'
#' @export
args = commandArgs(trailingOnly=TRUE)
intab <- read.delim(args[1], header=FALSE)
outtab <- reshape2::dcast(intab, V3 ~ factor(V1,levels=c("kingdom","phylum","class","order","family","genus","species")), value.var = "V2", na.rm = FALSE)
# errors in the previous command can occur if the levels of V1 are different
# solution 1
# need to see what the levels are:
unique(intab$V1)
# need to add the missing levels in the appropriate location on the hierarchy
outtab <- reshape2::dcast(intab, V3 ~ factor(V1,levels=c("kingdom","phylum","class","order","family","genus","species group","species","no rank")), value.var = "V2", na.rm = FALSE)

colnames(outtab)[1] <- "SpeciesID"
write.table(outtab, file = "wideformat.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#########################################################################

# 1. Replace all tabs with an underscore
tr '\t' '_' < all_taxa.txt > temp.txt

# 2. Replace all newlines with a double tab
tr '\n' '\t\t' < temp.txt > temp2.txt

# 3. Replace all occurrences of two consecutive tabs with a newline
sed 's/\t\t/\n/g' temp2.txt > output.txt

# 4 combine ID from blastN with output
awk -F, '{getline f1 <"blastN_out.tab" ;print f1,$0}' OFS=, output.txt > newout.txt


