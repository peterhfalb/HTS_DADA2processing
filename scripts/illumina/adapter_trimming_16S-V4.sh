#!/bin/bash
# =============================================================================
# Illumina 16S-V4 Adapter and Primer Trimming
# Removes TruSeq adapters and 515F/806R primers
# =============================================================================

set -euo pipefail

# Ensure modules are loaded
command -v cutadapt &> /dev/null || { echo "ERROR: cutadapt not loaded"; exit 1; }

echo "Removing adapters and primers from 16S-V4 Illumina data..."

# Use OUTPUT_DIR if set by SLURM script, otherwise use current directory
OUTPUT_DIR="${OUTPUT_DIR:-.}"
ADAPTER_DIR="$OUTPUT_DIR/01_adapter"
FILTERED_DIR="$OUTPUT_DIR/02_filtered"

# Create output directories
mkdir -p "$ADAPTER_DIR" "$FILTERED_DIR"

# Generate cutadapt commands for all samples
for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue  # Skip if R2 doesn't exist

  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 150 \
    -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA \
    -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT \
    -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT \
    -o $ADAPTER_DIR/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} \
    -p $ADAPTER_DIR/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} \
    ${i} ${i//_R1_/_R2_} \
    > $ADAPTER_DIR/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" \
  >> run_cutadapt.sh
done

chmod +x run_cutadapt.sh

echo "Running cutadapt trimming in parallel..."
parallel < run_cutadapt.sh

# Organize logs
mkdir -p "$ADAPTER_DIR/logs"
mv "$ADAPTER_DIR"/cutadapt*.log.txt "$ADAPTER_DIR/logs/" 2>/dev/null || true

# Summary
echo "Adapter trimming complete."
grep "passing" "$ADAPTER_DIR/logs/"* > "$OUTPUT_DIR/summary_adapter_trimming_16S-V4.txt" 2>/dev/null || echo "No logs found"

echo ""
echo "16S-V4 adapter trimming complete."
echo "Trimmed paired reads are in $ADAPTER_DIR/"
echo ""
