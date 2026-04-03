#!/bin/bash
# =============================================================================
# Illumina 18S-AMF (Arbuscular Mycorrhizal Fungi) Adapter and Primer Trimming
# Removes TruSeq adapters and WANDA_1/AML2_1 primers
# =============================================================================

set -euo pipefail

command -v cutadapt &> /dev/null || { echo "ERROR: cutadapt not loaded"; exit 1; }

echo "Removing adapters and primers from 18S-AMF Illumina data..."

mkdir -p 01_adapter

# WANDA_1 forward primer: CAGCCGCGGTAATTCCAGCT
# AML2_1 reverse primer: GAACCCAAACACTTTGGTTTCC

for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue

  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 \
    -g CAGCCGCGGTAATTCCAGCT \
    -G GAACCCAAACACTTTGGTTTCC \
    -o 01_adapter/${i//_R1_001.fastq.gz/_R1_001.fastq.gz} \
    -p 01_adapter/${i//_R1_001.fastq.gz/_R2_001.fastq.gz} \
    ${i} ${i//_R1_/_R2_} \
    > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" \
  >> run_cutadapt.sh
done

chmod +x run_cutadapt.sh
echo "Running cutadapt in parallel..."
parallel < run_cutadapt.sh

cd 01_adapter
mkdir -p logs
mv cutadapt*.log.txt logs/ 2>/dev/null || true

echo "18S-AMF adapter trimming complete."
grep "passing" logs/* > ../dada2output/summary_adapter_trimming_18S-AMF.txt 2>/dev/null || echo "Trimming complete"

cd ..
mkdir -p 02_filtered

echo ""
echo "18S-AMF adapter and primer trimming complete."
echo "Trimmed paired reads are in 01_adapter/"
echo ""
