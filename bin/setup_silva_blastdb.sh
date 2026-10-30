#!/bin/bash
set -euo pipefail

# One-time setup: build the SILVA SSU Ref NR99 BLAST database used by the
# 18S-AMF SILVA screen (rule amf_silva_filter).
# Usage: bash setup_silva_blastdb.sh [output_dir]
#
# Keep the default output path in sync with SILVA_BLAST_DB in the Snakefile
# and bin/run_dada2processing.

SILVA_VERSION="138.2"
SILVA_URL="https://www.arb-silva.de/fileadmin/silva_databases/release_${SILVA_VERSION}/Exports/SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta.gz"
OUT_DIR="${1:-/projects/standard/kennedyp/shared/taxonomy/silva_${SILVA_VERSION}_NR99_blastdb}"
DB_NAME="SILVA_${SILVA_VERSION}_NR99"
FASTA_GZ="$OUT_DIR/SILVA_${SILVA_VERSION}_SSURef_NR99_tax_silva.fasta.gz"
FASTA_DNA="$OUT_DIR/${DB_NAME}.fasta"

set +u
module load blast-plus/2.13.0-gcc-8.2.0-vo4mr4d 2>/dev/null || true
set -u
command -v makeblastdb &>/dev/null || { echo "ERROR: makeblastdb not found (module load blast-plus)"; exit 1; }

mkdir -p "$OUT_DIR"

if [[ ! -f "$FASTA_GZ" ]]; then
    echo "Downloading SILVA ${SILVA_VERSION} SSU Ref NR99..."
    wget -O "$FASTA_GZ" "$SILVA_URL"
fi

# SILVA sequences are RNA (U); convert to DNA on sequence lines only
echo "Converting U -> T..."
gunzip -c "$FASTA_GZ" | awk '/^>/ {print; next} {gsub(/U/, "T"); gsub(/u/, "t"); print}' > "$FASTA_DNA"

echo "Building BLAST database..."
makeblastdb -in "$FASTA_DNA" -dbtype nucl -title "$DB_NAME" -out "$OUT_DIR/$DB_NAME"

# Sanity check: AMF rank naming in this SILVA release
echo ""
echo "AMF taxonomy check (headers containing each name):"
for name in Glomeromycotina Glomeromycetes Glomeromycota Mucoromycota; do
    printf "  %-16s %s\n" "$name" "$(grep -c "^>.*$name" "$FASTA_DNA" || true)"
done
echo "Example AMF header:"
grep -m1 "^>.*Glomero" "$FASTA_DNA" || echo "  (none found)"

if ! grep -q "^>.*Glomeromycotina" "$FASTA_DNA"; then
    echo ""
    echo "WARNING: 'Glomeromycotina' does not appear in SILVA ${SILVA_VERSION} headers."
    echo "         Set amf_target_regex in the pipeline (Snakefile default) to a name listed above."
fi

rm -f "$FASTA_DNA"
echo ""
echo "Done. BLAST database: $OUT_DIR/$DB_NAME"
