#!/bin/bash

# remove adapters
module load cutadapt
module load parallel
module load trimmomatic/0.39

echo "creating output directories dada2output and 01_adapter"
mkdir dada2output
mkdir 01_adapter

echo "removing adapters"
for i in *_R1_001.fastq.gz; do echo "java -jar /users/4/goul0109/Trimmomatic-0.39/trimmomatic-0.39.jar PE $i ${i//_R1_/_R2_} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} ILLUMINACLIP:/users/4/goul0109/scripts/aviti_adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100 -trimlog 01_adapter/${i//_R1_001.fastq.gz/_log.txt}" >> run_trim.sh; done
chmod +x run_trim.sh
parallel < run_trim.sh
cd 01_adapter
mkdir 01_logs
mv *_log.txt 01_logs
rm *unpaired.fastq.gz
grep "passing" 01_logs/* > ../dada2output/summary_adapter_trimming.txt
echo "adapters removed"

echo "removing V4 primers"
# remove primers
mkdir ../02_filtered
for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -a GTGCCAGCMGCCGCGGTAA...ATTAGANACCCNNGTAGTCC -A GGACTACHVGGGTWTCTAAT...TTACCGCGGCNGCTGGCAC -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > ../dada2output/summary_primer_trimming.txt
echo "primers removed"

echo "removing extra 16S primers"
mkdir ../03_filtered
for i in *_R1_001.fastq.gz; do echo "java -jar /users/4/goul0109/Trimmomatic-0.39/trimmomatic-0.39.jar PE $i ${i//_R1/_R2} ../03_filtered/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R1.unpaired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R2.unpaired.fastq.gz} ILLUMINACLIP:/users/4/goul0109/scripts/aviti_primer_as_adapter.fa:2:30:10:2:True MINLEN:100 -trimlog ../03_filtered/${i//_R1_001.fastq.gz/_log.txt}" >> re_adapt.sh; done
chmod +x re_adapt.sh
parallel < re_adapt.sh
cd ../03_filtered/
rm *unpaired.fastq.gz
mkdir 03_logs
mv *log.txt 03_logs
grep "passing" 03_logs/* > ../dada2output/summary_primer_trimming.txt
echo "extra primers removed"

