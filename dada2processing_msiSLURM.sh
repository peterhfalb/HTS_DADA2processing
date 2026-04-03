#!/bin/bash
# =============================================================================
# DADA2 Processing Pipeline - SLURM Job Script (MSI AGATE)
# Submitted via run_dada2processing wrapper
# =============================================================================

#SBATCH --job-name=DADA2pipeline
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=60gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err

set -euo pipefail

# Load user configuration (PIPELINE_DIR passed by wrapper)
[ -n "${PIPELINE_DIR:-}" ] || { echo "ERROR: PIPELINE_DIR not set. Submit via run_dada2processing wrapper."; exit 1; }

# Load required modules
module purge
module load R/4.4.0
module load cutadapt
module load trimmomatic

echo "✓ Required modules loaded"
echo "R: $(Rscript --version 2>&1)"

# Parse arguments
if [ "$#" -ne 5 ]; then
  echo "Usage: dada2processing_msiSLURM.sh <input_dir> <output_dir> <marker_gene> <platform> <quality>"
  exit 1
fi

INPUT_DIR="$1"
OUTPUT_DIR="$2"
MARKER_GENE="$3"
PLATFORM="$4"
QUALITY="$5"

# Setup output directories
mkdir -p "$OUTPUT_DIR/logs"
cd "$INPUT_DIR"

echo "=========================================="
echo "DADA2 Processing Pipeline"
echo "=========================================="
echo "Start time: $(date)"
echo "Input directory:  $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Marker gene:      $MARKER_GENE"
echo "Platform:         $PLATFORM"
echo "Quality:          $QUALITY"
echo "=========================================="
echo ""

# Define expected primers for validation
declare -A EXPECTED_PRIMERS
EXPECTED_PRIMERS["16S-V4_F"]="GTGYCAGCMGCCGCGGTAA"
EXPECTED_PRIMERS["16S-V4_R"]="GGACTACNVGGGTWTCTAAT"
EXPECTED_PRIMERS["ITS1_F"]="CTTGGTCATTTAGAGGAAGTAA"
EXPECTED_PRIMERS["ITS1_R"]="GCTGCGTTCTTCATCGATGC"
EXPECTED_PRIMERS["ITS2_F"]="TCGATGAAGAACGCAGCG"
EXPECTED_PRIMERS["ITS2_R"]="TCCTCCGCTTATTGATATGC"
EXPECTED_PRIMERS["18S-V4_F"]="TTAAARVGYTCGTAGTYG"
EXPECTED_PRIMERS["18S-V4_R"]="CCGTCAATTHCTTYAART"
EXPECTED_PRIMERS["18S-AMF_F"]="CAGCCGCGGTAATTCCAGCT"
EXPECTED_PRIMERS["18S-AMF_R"]="GAACCCAAACACTTTGGTTTCC"

# Step 1: Primer validation
echo ""
echo "Step 1: Validating primers..."
echo "========================================"

FIRST_R1=$(ls *_R1*.fastq.gz 2>/dev/null | head -1)
FIRST_R2=$(ls *_R2*.fastq.gz 2>/dev/null | head -1)

if [ -z "$FIRST_R1" ] || [ -z "$FIRST_R2" ]; then
  echo "ERROR: Could not find paired FASTQ files in $INPUT_DIR"
  exit 1
fi

# Run primercheck (stub for now — detailed validation will happen in actual primercheck script)
# This is a sanity check to detect gross primer mismatches
"$PIPELINE_DIR/scripts/primercheck.sh" "$FIRST_R1" "$FIRST_R2" | tee "$OUTPUT_DIR/logs/primercheck.log"

echo "Primer check complete."
echo ""

# Step 2: Adapter trimming
echo "Step 2: Adapter and primer trimming..."
echo "========================================"

# Export variables for adapter trimming scripts
export OUTPUT_DIR
export INPUT_DIR

case "$MARKER_GENE,$PLATFORM" in
  16S-V4,illumina)
    bash "$PIPELINE_DIR/scripts/illumina/adapter_trimming_16S-V4.sh"
    ;;
  16S-V4,aviti)
    bash "$PIPELINE_DIR/scripts/aviti/adapter_trimming_16S-V4.sh"
    ;;
  ITS1,illumina)
    bash "$PIPELINE_DIR/scripts/illumina/adapter_trimming_ITS1.sh"
    ;;
  ITS1,aviti)
    bash "$PIPELINE_DIR/scripts/aviti/adapter_trimming_ITS1.sh"
    ;;
  ITS2,illumina)
    bash "$PIPELINE_DIR/scripts/illumina/adapter_trimming_ITS2.sh"
    ;;
  ITS2,aviti)
    bash "$PIPELINE_DIR/scripts/aviti/adapter_trimming_ITS2.sh"
    ;;
  18S-V4,illumina)
    bash "$PIPELINE_DIR/scripts/illumina/adapter_trimming_18S-V4.sh"
    ;;
  18S-V4,aviti)
    bash "$PIPELINE_DIR/scripts/aviti/adapter_trimming_18S-V4.sh"
    ;;
  18S-AMF,illumina)
    bash "$PIPELINE_DIR/scripts/illumina/adapter_trimming_18S-AMF.sh"
    ;;
  18S-AMF,aviti)
    bash "$PIPELINE_DIR/scripts/aviti/adapter_trimming_18S-AMF.sh"
    ;;
  *)
    echo "ERROR: Unknown marker_gene/platform combination: $MARKER_GENE/$PLATFORM"
    exit 1
    ;;
esac

echo "Adapter and primer trimming complete."
echo ""

# Step 3: DADA2 processing
echo "Step 3: DADA2 ASV inference..."
echo "========================================"

# The filtered reads should be in 02_filtered/ (after adapter trimming)
if [ ! -d "02_filtered" ]; then
  echo "ERROR: 02_filtered directory not found. Adapter trimming may have failed."
  exit 1
fi

cd 02_filtered

case "$MARKER_GENE,$PLATFORM" in
  16S-V4,illumina)
    Rscript "$PIPELINE_DIR/scripts/illumina/run_16S-V4_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  16S-V4,aviti)
    Rscript "$PIPELINE_DIR/scripts/aviti/run_16S-V4_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  ITS1,illumina)
    Rscript "$PIPELINE_DIR/scripts/illumina/run_ITS_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  ITS1,aviti)
    Rscript "$PIPELINE_DIR/scripts/aviti/run_ITS_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  ITS2,illumina)
    Rscript "$PIPELINE_DIR/scripts/illumina/run_ITS_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  ITS2,aviti)
    Rscript "$PIPELINE_DIR/scripts/aviti/run_ITS_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  18S-V4,illumina)
    Rscript "$PIPELINE_DIR/scripts/illumina/run_18S-V4_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  18S-V4,aviti)
    Rscript "$PIPELINE_DIR/scripts/aviti/run_18S-V4_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  18S-AMF,illumina)
    Rscript "$PIPELINE_DIR/scripts/illumina/run_18S-AMF_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  18S-AMF,aviti)
    Rscript "$PIPELINE_DIR/scripts/aviti/run_18S-AMF_dada2.R" "$QUALITY" "$OUTPUT_DIR"
    ;;
  *)
    echo "ERROR: Unknown marker_gene/platform combination: $MARKER_GENE/$PLATFORM"
    exit 1
    ;;
esac

echo "DADA2 processing complete."
echo ""

# Write summary
echo "=========================================="
echo "Pipeline Completion Summary"
echo "=========================================="
echo "End time: $(date)"
echo "Marker gene: $MARKER_GENE"
echo "Platform: $PLATFORM"
echo "Quality: $QUALITY"
echo ""
echo "Output directory: $OUTPUT_DIR"
echo "Check logs: ls -lh $OUTPUT_DIR/logs/"
echo "=========================================="

echo "PIPELINE COMPLETED SUCCESSFULLY" > "$OUTPUT_DIR/pipeline_run.log"
