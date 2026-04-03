#!/bin/bash

# remove adapters
module load cutadapt
module load parallel

echo "creating output directory 01_adapter"
mkdir 01_adapter

for i in *_R1_001.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 50 -n 5 -a TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -G GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -o 01_adapter/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh

echo "running cutadapt trimming nextera transposase adapters "
parallel < run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv cutadapt* 01_logs

grep "passing" 01_logs/* > ../dada2output/summary_adapter_trimming.txt
echo "summary_adapter_trimming.txt created"

echo "adapters removed"
