#!/bin/bash
# =============================================================================
# Illumina ITS1 Adapter and Primer Trimming
# Removes TruSeq adapters and ITS1F/ITS2 primers
# =============================================================================

set -euo pipefail

command -v cutadapt &> /dev/null || { echo "ERROR: cutadapt not loaded"; exit 1; }

echo "Removing adapters and primers from ITS1 Illumina data..."

mkdir -p 01_adapter

# ITS1F forward primer: CTTGGTCATTTAGAGGAAGTAA
# ITS2 reverse primer: GCTGCGTTCTTCATCGATGC

for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue

  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 \
    -g CTTGGTCATTTAGAGGAAGTAA \
    -G GCTGCGTTCTTCATCGATGC \
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

echo "ITS1 adapter trimming complete."
grep "passing" logs/* > ../dada2output/summary_adapter_trimming_ITS1.txt 2>/dev/null || echo "Trimming complete"

cd ..
mkdir -p 02_filtered

echo ""
echo "ITS1 adapter and primer trimming complete."
echo "Trimmed paired reads are in 01_adapter/"
echo ""
