#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=64
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=goul0109@umn.edu

cd /location/of/fastq.gz/files/

module load parallel
module load cutadapt
module load R/4.4.0-openblas-rocky8

mkdir dada2output
adapter_removal.sh
# what are the primers
primer_trim.sh -p GTGCCAGCMGCCGCGGTAA -q GGACTACHVGGGTWTCTAAT
# what is the quality
run_dada2.sh 16S good