#!/bin/bash
# =============================================================================
# DADA2 Primer Detection and Validation Script
# Checks which primers are present in FASTQ files
# Validates that they match the specified marker gene
# =============================================================================

set -euo pipefail

print_help() {
cat << 'HELP'
primercheck.sh — Detect and validate primers in FASTQ files

USAGE:
  primercheck.sh <R1_fastq> <R2_fastq> [<expected_marker_gene>]

ARGUMENTS:
  R1_fastq              Path to forward (R1) FASTQ file (can be .gz)
  R2_fastq              Path to reverse (R2) FASTQ file (can be .gz)
  expected_marker_gene  Optional: expected marker gene for validation
                        16S-V4, ITS1, ITS2, 18S-V4, 18S-AMF

OUTPUT:
  Detected marker genes and primers
  Validation result (OK or WARNING)

HELP
}

if [ $# -lt 2 ] || [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
  print_help
  exit 0
fi

R1_FILE="$1"
R2_FILE="$2"
EXPECTED_MARKER="${3:-}"

if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
  echo "ERROR: FASTQ files not found"
  exit 1
fi

# Define UMGC primer sequences (exact matches, not regex)
declare -A PRIMER_R1
PRIMER_R1["16S-V4"]="GTGYCAGCMGCCGCGGTAA"
PRIMER_R1["ITS1"]="CTTGGTCATTTAGAGGAAGTAA"
PRIMER_R1["ITS2"]="TCGATGAAGAACGCAGCG"
PRIMER_R1["18S-V4"]="TTAAARVGYTCGTAGTYG"
PRIMER_R1["18S-AMF"]="CAGCCGCGGTAATTCCAGCT"

declare -A PRIMER_R2
PRIMER_R2["16S-V4"]="GGACTACNVGGGTWTCTAAT"
PRIMER_R2["ITS1"]="GCTGCGTTCTTCATCGATGC"
PRIMER_R2["ITS2"]="TCCTCCGCTTATTGATATGC"
PRIMER_R2["18S-V4"]="CCGTCAATTHCTTYAART"
PRIMER_R2["18S-AMF"]="GAACCCAAACACTTTGGTTTCC"

# Function to search for primers (with degenerate bases)
search_primer() {
  local file="$1"
  local primer="$2"
  local max_count=50  # Search first 50 matches

  # Convert degenerate bases to regex character class for grep
  # This is a simplified conversion for common degenerate codes
  local grep_pattern="$primer"
  grep_pattern="${grep_pattern//Y/[CT]}"
  grep_pattern="${grep_pattern//R/[AG]}"
  grep_pattern="${grep_pattern//W/[AT]}"
  grep_pattern="${grep_pattern//S/[GC]}"
  grep_pattern="${grep_pattern//K/[GT]}"
  grep_pattern="${grep_pattern//M/[AC]}"
  grep_pattern="${grep_pattern//B/[CGT]}"
  grep_pattern="${grep_pattern//D/[AGT]}"
  grep_pattern="${grep_pattern//H/[ACT]}"
  grep_pattern="${grep_pattern//V/[ACG]}"
  grep_pattern="${grep_pattern//N/.}"

  # Search in the FASTQ file (handle gzip and uncompressed)
  local count=0
  if [[ "$file" == *.gz ]]; then
    count=$(zgrep -c "$grep_pattern" "$file" 2>/dev/null || echo 0)
  else
    count=$(grep -c "$grep_pattern" "$file" 2>/dev/null || echo 0)
  fi

  echo "$count"
}

echo "=========================================="
echo "Primer Detection"
echo "=========================================="
echo "R1 file: $R1_FILE"
echo "R2 file: $R2_FILE"
echo ""

DETECTED_MARKERS=""
DETECTION_THRESHOLD=10  # At least 10 reads should have the primer

for marker in "16S-V4" "ITS1" "ITS2" "18S-V4" "18S-AMF"; do
  R1_PRIMER="${PRIMER_R1[$marker]}"
  R2_PRIMER="${PRIMER_R2[$marker]}"

  R1_COUNT=$(search_primer "$R1_FILE" "$R1_PRIMER")
  R2_COUNT=$(search_primer "$R2_FILE" "$R2_PRIMER")

  if [ "$R1_COUNT" -ge "$DETECTION_THRESHOLD" ] && [ "$R2_COUNT" -ge "$DETECTION_THRESHOLD" ]; then
    echo "✓ Detected: $marker (R1: $R1_COUNT, R2: $R2_COUNT)"
    DETECTED_MARKERS="$DETECTED_MARKERS $marker"
  fi
done

echo ""

if [ -z "$DETECTED_MARKERS" ]; then
  echo "WARNING: No standard UMGC primers detected in the FASTQ files."
  echo "This may indicate:"
  echo "  1. Primers have been removed already"
  echo "  2. The files are empty or corrupted"
  echo "  3. The data uses non-standard primers"
  echo ""
fi

# Validation against expected marker
if [ -n "$EXPECTED_MARKER" ]; then
  echo "Validation: Expected marker = $EXPECTED_MARKER"
  if [[ "$DETECTED_MARKERS" == *"$EXPECTED_MARKER"* ]]; then
    echo "✓ PASS: Detected marker matches expected marker"
  else
    echo "⚠ WARNING: Detected markers ($DETECTED_MARKERS) do not include expected marker ($EXPECTED_MARKER)"
    echo "    This may cause processing issues. Please verify the marker gene specification."
  fi
fi

echo "=========================================="
