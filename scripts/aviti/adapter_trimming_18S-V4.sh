#!/bin/bash
# Aviti 18S-V4 (Microeukaryote) Adapter and Primer Trimming (3-step)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AVITI_ADAPTERS="$SCRIPT_DIR/aviti_adapters.fa"
AVITI_PRIMERS="$SCRIPT_DIR/aviti_ITSprimer_as_adapter.fa"

[ -f "$AVITI_ADAPTERS" ] && [ -f "$AVITI_PRIMERS" ] || { echo "ERROR: Aviti adapter FASTA files not found"; exit 1; }

mkdir -p 01_adapter 02_filtered 03_filtered

echo "18S-V4 Aviti trimming..."
for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue
  echo "trimmomatic PE $i ${i//_R1_/_R2_} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} 01_adapter/${i//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} ILLUMINACLIP:$AVITI_ADAPTERS:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:100 -trimlog 01_adapter/${i//_R1_001.fastq.gz/_log.txt}" >> run_trim.sh
done
chmod +x run_trim.sh && parallel < run_trim.sh
cd 01_adapter && mkdir -p logs && mv *_log.txt logs/ 2>/dev/null || true && rm *unpaired.fastq.gz 2>/dev/null || true && cd ..

for i in 01_adapter/*_R1_001.paired.fastq.gz; do
  basename="${i##*/}"
  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -g TTAAARVGYTCGTAGTYG -G CCGTCAATTHCTTYAART -o 02_filtered/${basename//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p 02_filtered/${basename//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 02_filtered/cutadapt.${basename//_R1_001.paired.fastq.gz/.log.txt}" >> run_cutadapt.sh
done
chmod +x run_cutadapt.sh && parallel < run_cutadapt.sh

for i in 02_filtered/*_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue
  basename="${i##*/}"
  echo "trimmomatic PE ${i} ${i//_R1_/_R2_} 03_filtered/${basename//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} 03_filtered/${basename//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} 03_filtered/${basename//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} 03_filtered/${basename//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} ILLUMINACLIP:$AVITI_PRIMERS:2:30:10:2:True MINLEN:100 -trimlog 03_filtered/${basename//_R1_001.fastq.gz/_log.txt}" >> run_trim2.sh
done
chmod +x run_trim2.sh && parallel < run_trim2.sh
cd 03_filtered && rm *unpaired.fastq.gz 2>/dev/null || true && cp *_R1_001.paired.fastq.gz ../02_filtered/ 2>/dev/null || true && cp *_R2_001.paired.fastq.gz ../02_filtered/ 2>/dev/null || true && cd ../02_filtered
for f in *_R1_001.paired.fastq.gz; do mv "$f" "${f//_R1_001.paired.fastq.gz/_R1_001.fastq.gz}" 2>/dev/null || true; done
for f in *_R2_001.paired.fastq.gz; do mv "$f" "${f//_R2_001.paired.fastq.gz/_R2_001.fastq.gz}" 2>/dev/null || true; done
echo "18S-V4 Aviti trimming complete."
