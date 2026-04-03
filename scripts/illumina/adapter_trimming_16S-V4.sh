#!/bin/bash
# =============================================================================
# Illumina 16S-V4 Adapter and Primer Trimming
# Removes TruSeq adapters and 515F/806R primers
# =============================================================================

set -euo pipefail

# Ensure modules are loaded
command -v cutadapt &> /dev/null || { echo "ERROR: cutadapt not loaded"; exit 1; }

echo "Removing adapters and primers from 16S-V4 Illumina data..."

# Create output directory
mkdir -p 01_adapter

# Generate cutadapt commands for all samples
for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue  # Skip if R2 doesn't exist

  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 \
    -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA \
    -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT \
    -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT \
    -o 01_adapter/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} \
    -p 01_adapter/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} \
    ${i} ${i//_R1_/_R2_} \
    > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" \
  >> run_cutadapt.sh
done

chmod +x run_cutadapt.sh

echo "Running cutadapt trimming in parallel..."
parallel < run_cutadapt.sh

# Organize logs
cd 01_adapter
mkdir -p logs
mv cutadapt*.log.txt logs/ 2>/dev/null || true

# Summary
echo "Adapter trimming complete."
grep "passing" logs/* > ../dada2output/summary_adapter_trimming_16S-V4.txt 2>/dev/null || echo "No logs found"

# Prepare for primer trimming (next step)
echo "Creating 02_filtered directory for primer trimming..."
cd ..
mkdir -p 02_filtered

# Note: Primers are trimmed in the above cutadapt call using the -g and -a flags
# 515F: GTGYCAGCMGCCGCGGTAA (forward, at 5' of R1)
# 806R: GGACTACNVGGGTWTCTAAT (reverse, at 5' of R2 after reverse complement)
# These will be trimmed in the next step via cutadapt with primer-specific command

echo ""
echo "16S-V4 adapter trimming complete."
echo "Trimmed paired reads are in 01_adapter/"
echo ""
