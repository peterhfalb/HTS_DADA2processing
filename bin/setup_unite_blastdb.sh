#!/bin/bash
set -euo pipefail

# One-time setup: build a BLAST database from the lab's UNITE fasta, used by
# check_primers --blast (without it, check_primers builds a temporary one each run).
# Usage: bash setup_unite_blastdb.sh [output_dir]
#
# Keep paths in sync with UNITE_FASTA / UNITE_DB in bin/check_primers.

TAX_DIR="/projects/standard/kennedyp/shared/taxonomy"
UNITE_FASTA="$TAX_DIR/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta"
OUT_DIR="${1:-$TAX_DIR/unite_blastdb}"
DB_NAME="UNITE_all_19.02.2025"

set +u
module load blast-plus/2.13.0-gcc-8.2.0-vo4mr4d 2>/dev/null || true
set -u
command -v makeblastdb &>/dev/null || { echo "ERROR: makeblastdb not found (module load blast-plus)"; exit 1; }
[[ -f "$UNITE_FASTA" ]] || { echo "ERROR: UNITE fasta not found at $UNITE_FASTA"; exit 1; }

mkdir -p "$OUT_DIR"
echo "Building BLAST database..."
makeblastdb -in "$UNITE_FASTA" -dbtype nucl -title "$DB_NAME" -out "$OUT_DIR/$DB_NAME"

echo ""
echo "Done. BLAST database: $OUT_DIR/$DB_NAME"
