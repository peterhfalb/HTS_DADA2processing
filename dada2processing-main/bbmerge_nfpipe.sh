#!/bin/bash -l
#SBATCH --time=48:00:00
#SBATCH --ntasks=128
#SBATCH --mem=200gb
#SBATCH --mail-type=ALL
#SBATCH --account=moeller
#SBATCH --mail-user=goul0109@umn.edu

cd /scratch.global/goul0109/moeller2/03_filtered/filtered/
mkdir ../nf_pipe

for i in *_F_filt.fastq.gz; do mkdir ../nf_pipe/${i//_F_filt.fastq.gz/}: done
for i in *_F_filt.fastq.gz; do echo "/home/umii/goul0109/bbmap/bbmerge.sh in1=$i in2=${i//_F/_R} out=../nf_pipe/${i//_F_filt.fastq.gz/}/${i//_F_filt.fastq.gz/.fastq.gz}" >> run_bbmerge.cmd; done
chmod +x run_bbmerge.cmd
./run_bbmerge.cmd

module load java/openjdk-21.0.2
module load singularity
cd /scratch.global/goul0109/moeller2/03_filtered/
nextflow run epi2me-labs/wf-metagenomics --fastq nf_pipe --classifier kraken2 --database_set ncbi_16s_18s -profile singularity