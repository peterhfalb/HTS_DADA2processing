#!/bin/bash
# usage
# To call the script:
# bash primer_trim.sh -p forwardprimer -q reverseprimer
# primer_trim.sh -p GTGCCAGCMGCCGCGGTAA -q GGACTACHVGGGTWTCTAAT
while getopts "p:q:" OPTION;
do
        case $OPTION in
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
mkdir ../02_filtered  
for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g ${forward} -G ${reverse} -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > ../dada2output/summary_primer_trimming.txt
echo "primers trimmed"
