#!/bin/bash
# single end fasta trim intended for taxonomy database trimming
# requires bbtools in your PATH
# usage
# To call the script:
# bash fasta_primer_trim.sh -p forwardprimer -q reverseprimer
# primer_trim.sh -f file.fasta -p GTGCCAGCMGCCGCGGTAA -q GGACTACHVGGGTWTCTAAT
while getopts "f:p:q:" OPTION;
do
        case $OPTION in
        f)
        fasta=$OPTARG
        ;;
        p)
        forward=$OPTARG
        ;;
        q)
        reverse=$OPTARG
        ;;
  esac
done
#primer trimming
echo "Primer 1: ${forward}"
echo "Primer 2: ${reverse}"

# remove adapters
module load cutadapt
module load parallel

# remove primers
cutadapt --cores 8 --minimum-length 100 --discard-untrimmed -g ${forward} -G ${reverse} -o exact.fasta ${fasta} > cutadapt.primer.log.txt

#needs to have this format forward_primer...reverse_comp_rev_primer
cutadapt --cores 8 --minimum-length 100 --discard-untrimmed -g GTGYCAGCMGCCGCGGTAA...ATTAGANACCCNNGTAGTCC -o exact.fasta /home/umii/public/dada2_taxonomy_references/silva_nr99_v138.1_train_set.fa > cutadapt.primer.log.txt

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd

mkdir primer_logs
mv *log.txt primer_logs
grep "passing" primer_logs/* > summary_primer_trimming.txt
echo "primers trimmed from fasta"

# dedupe

/home/umii/goul0109/bbmap/dedupe.sh in=exact.fasta out=dedupe.fasta mergenames=t exact=t
#/home/umii/goul0109/bbmap/

echo "sequences deduplicated and fasta headers merged"

# headers will now look like this with > separating taxonomy of merged sequences
# >Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Saccharopolyspora;>Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Pseudonocardia;
# to
# >Bacteria;Actinobacteria;Actinobacteria;Pseudonocardiales;Pseudonocardiaceae;Saccharopolyspora;>;Pseudonocardia;

awk -F '\t' '{delete seen; c=0; for (i=1;i<=NF;i++) if (!seen[$i]++) printf "%s%s", (++c>1?OFS:""), $i; print ""}' dedupe.fasta > exactDB.fasta
# combine separate headers
awk  -F '>' '{OFS = FS} {delete seen; c=0; for (i=1;i<=NF;i++) if (!seen[$i]++) printf "%s%s", (++c>1?OFS:""), $i; print ""}' exactDB.fasta > tempexactDB1.fasta
# combine taxa levels
awk  -F ';' '{OFS = FS} {delete seen; c=0; for (i=1;i<=NF;i++) if (!seen[$i]++) printf "%s%s", (++c>1?OFS:""), $i; print ""}' tempexactDB1.fasta > tempexactDB2.fasta
#combine different taxa at genus level with _
sed -i ':a;/[;].*[;]/s/[;]/_/6;ta' tempexactDB2.fasta
sed -i -E 's/(.*),/\1;/' tempexactDB2.fasta
# if less than 5 columns exist it doesn't keep the ; at the end, which is bad. 
sed -i '/;/ s/$/;/' tempexactDB2.fasta
sed -e 's/;;/;/g' tempexactDB2.fasta > FinalExactDB.fasta
# complete message
echo "output complete: exactDB_${fasta}"
