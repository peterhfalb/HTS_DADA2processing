#!/bin/bash

REGION=$1

echo "creating output directories dada2output and 01_adapter"
mkdir dada2output
mkdir 01_adapter

echo "removing adapters"
for i in *_R1_001.fastq.gz; do echo "java -jar /home/umii/goul0109/Trimmomatic-0.39/trimmomatic-0.39.jar PE $i ${i//_R1_/_R2_} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} ILLUMINACLIP:/home/umii/goul0109/scripts/aviti_adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100 -trimlog 01_adapter/${i//_R1_001.fastq.gz/_log.txt}" >> run_trim.sh; done
chmod +x run_trim.sh
parallel < run_trim.sh
cd 01_adapter    
mkdir 01_logs
mv *_log.txt 01_logs
rm *unpaired.fastq.gz
echo "adapters removed"

echo "removing ITS primers"
# remove primers
mkdir ../02_filtered

# REGION must be one of: ITS1, ITS2, V4, V3V5, V5V6, V3V4

/scratch.global/goul0109/dada2processing/generate_cutadapt.sh $REGION

chmod +x run_cutadapt2.cmd
parallel < run_cutadapt2.cmd
cd ../02_filtered/
mkdir 02_logs
mv *log.txt 02_logs
grep "passing" 02_logs/* > ../dada2output/summary_primer_trimming.txt
echo "primers removed"

echo "removing extra ITS primers"
mkdir ../03_filtered
for i in *_R1_001.fastq.gz; do echo "java -jar /home/umii/goul0109/Trimmomatic-0.39/trimmomatic-0.39.jar PE $i ${i//_R1/_R2} ../03_filtered/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R1.unpaired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} ../03_filtered/${i//_R1_001.fastq.gz/_R2.unpaired.fastq.gz} ILLUMINACLIP:/home/umii/goul0109/scripts/aviti_ITSprimer_as_adapter.fa:2:30:10:2:True MINLEN:100 -trimlog ../03_filtered/${i//_R1_001.fastq.gz/_log.txt}" >> re_adapt.sh; done
chmod +x re_adapt.sh
parallel < re_adapt.sh
cd ../03_filtered/
rm *unpaired.fastq.gz
mkdir 03_logs
mv *log.txt 03_logs
echo "extra primers removed"


