# dada2processing
 process from a directory of fastq.gz files to a dada2 table
 
 ## Install
 git clone https://github.com/trevorjgould/dada2processing.git
 
 add dada2processing directory to your $PATH
 
 ## create output directory
 mkdir dada2output
 
 ## removes adapters
 adapter_removal.sh
 
 ## removes primers
 checks what the primers are: give an example pair of fastq.gz files
 
 primercheck.sh sample_R1_L001.fastq.gz sample_R2_L001.fastq.gz
 
 automatically runs: primer_trim.sh which removes the found primers
 	
 ## Run dada2: what markergene?, what is the quality?
 run_dada2.sh [-m < 16S|ITS|18S>] [-q <good|bad>]
 
 runs: run_16S_dada2.R / run_ITS_dada2.R / run_18S_dada2.R

## Here is an example slurm script that does all of this:

```
#!/bin/bash -l        
#SBATCH --time=12:00:00
#SBATCH --ntasks=64
#SBATCH --mem=60gb
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=USER@umn.edu

cd /location/of/fastq.gz/files/

module load parallel
module load cutadapt
module load R/4.4.0-openblas-rocky8

mkdir dada2output

adapter_removal.sh
cd 01_adapter

# what are the primers
primercheck.sh sample_R1_L001.fastq.gz sample_R2_L001.fastq.gz

cd 02_filtered

# what is the quality
run_dada2.sh -m 16S -q good
```