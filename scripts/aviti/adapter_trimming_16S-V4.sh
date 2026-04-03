#!/bin/bash
# =============================================================================
# Aviti 16S-V4 Adapter and Primer Trimming (3-step process)
# Step 1: Trimmomatic - removes Aviti adapters
# Step 2: cutadapt - removes 515f/806r primers
# Step 3: Trimmomatic - catches read-through primers
# =============================================================================

set -euo pipefail

# Check for required adapter FASTA files
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AVITI_ADAPTERS="$SCRIPT_DIR/aviti_adapters.fa"
AVITI_PRIMERS="$SCRIPT_DIR/aviti_ITSprimer_as_adapter.fa"

if [ ! -f "$AVITI_ADAPTERS" ] || [ ! -f "$AVITI_PRIMERS" ]; then
  echo "ERROR: Aviti adapter FASTA files not found in $SCRIPT_DIR"
  echo "  Missing: $([ ! -f "$AVITI_ADAPTERS" ] && echo "aviti_adapters.fa") $([ ! -f "$AVITI_PRIMERS" ] && echo "aviti_ITSprimer_as_adapter.fa")"
  echo "Please request these files from Trevor Gould."
  exit 1
fi

echo "Removing adapters and primers from 16S-V4 Aviti data..."
echo "Using adapter files from: $SCRIPT_DIR"

mkdir -p 01_adapter 02_filtered 03_filtered

# Step 1: Remove Aviti platform adapters with Trimmomatic
echo "Step 1: Removing Aviti adapters..."
for i in *_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue
  echo "java -jar \$TRIMMOMATIC_JAR PE $i ${i//_R1_/_R2_} \
    01_adapter/${i//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} \
    01_adapter/${i//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} \
    01_adapter/${i//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} \
    01_adapter/${i//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} \
    ILLUMINACLIP:$AVITI_ADAPTERS:2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:100 \
    -trimlog 01_adapter/${i//_R1_001.fastq.gz/_log.txt}" >> run_trim.sh
done

chmod +x run_trim.sh
TRIMMOMATIC_JAR=$(find /panfs -name "trimmomatic-0.39.jar" 2>/dev/null | head -1)
if [ -z "$TRIMMOMATIC_JAR" ]; then
  TRIMMOMATIC_JAR="trimmomatic-0.39.jar"  # Assume in PATH
fi
sed -i "s|\\\$TRIMMOMATIC_JAR|$TRIMMOMATIC_JAR|g" run_trim.sh
parallel < run_trim.sh

cd 01_adapter
mkdir -p logs
mv *_log.txt logs/ 2>/dev/null || true
rm *unpaired.fastq.gz 2>/dev/null || true

echo "Step 1 complete. Adapters removed."

# Step 2: Remove 515f/806r primers with cutadapt
echo "Step 2: Removing 515f/806r primers..."
cd ..
for i in 01_adapter/*_R1_001.paired.fastq.gz; do
  basename="${i##*/}"
  echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed \
    -g GTGYCAGCMGCCGCGGTAA \
    -G GGACTACNVGGGTWTCTAAT \
    -o 02_filtered/${basename//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} \
    -p 02_filtered/${basename//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} \
    ${i} ${i//_R1_/_R2_} \
    > 02_filtered/cutadapt.${basename//_R1_001.paired.fastq.gz/.log.txt}" >> run_cutadapt.sh
done

chmod +x run_cutadapt.sh
parallel < run_cutadapt.sh

cd 02_filtered
mkdir -p logs
mv *.log.txt logs/ 2>/dev/null || true
echo "Step 2 complete. Primers removed."

# Step 3: Remove read-through primers with Trimmomatic (optional, for extra safety)
echo "Step 3: Catching read-through primers with Trimmomatic..."
cd ..
for i in 02_filtered/*_R1_001.fastq.gz; do
  [ -f "${i//_R1_/_R2_}" ] || continue
  basename="${i##*/}"
  echo "java -jar $TRIMMOMATIC_JAR PE ${i} ${i//_R1_/_R2_} \
    03_filtered/${basename//_R1_001.fastq.gz/_R1_001.paired.fastq.gz} \
    03_filtered/${basename//_R1_001.fastq.gz/_R1_001.unpaired.fastq.gz} \
    03_filtered/${basename//_R1_001.fastq.gz/_R2_001.paired.fastq.gz} \
    03_filtered/${basename//_R1_001.fastq.gz/_R2_001.unpaired.fastq.gz} \
    ILLUMINACLIP:$AVITI_PRIMERS:2:30:10:2:True MINLEN:100 \
    -trimlog 03_filtered/${basename//_R1_001.fastq.gz/_log.txt}" >> run_trim2.sh
done

chmod +x run_trim2.sh
parallel < run_trim2.sh

cd 03_filtered
rm *unpaired.fastq.gz 2>/dev/null || true
mkdir -p logs
mv *_log.txt logs/ 2>/dev/null || true

# Copy final filtered reads back to 02_filtered for DADA2 processing
cp *_R1_001.paired.fastq.gz ../02_filtered/ 2>/dev/null || true
cp *_R2_001.paired.fastq.gz ../02_filtered/ 2>/dev/null || true

cd ../02_filtered
# Rename paired files to match DADA2 expectations
for f in *_R1_001.paired.fastq.gz; do
  mv "$f" "${f//_R1_001.paired.fastq.gz/_R1_001.fastq.gz}" 2>/dev/null || true
done
for f in *_R2_001.paired.fastq.gz; do
  mv "$f" "${f//_R2_001.paired.fastq.gz/_R2_001.fastq.gz}" 2>/dev/null || true
done

echo ""
echo "16S-V4 Aviti adapter trimming complete (3-step process)."
echo "Final filtered reads are in 02_filtered/"
echo ""
