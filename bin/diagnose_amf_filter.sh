#!/bin/bash
set -euo pipefail

# Diagnostics for the 18S-AMF SILVA screen (rule amf_silva_filter).
# Breaks down discarded OTUs/reads by reason, then BLASTs the most abundant discarded
# OTUs against SILVA, MaarjAM, Eukaryome (and optionally ROD) without the identity and
# coverage cutoffs, to show why they failed the screen.
#
# Usage: bash diagnose_amf_filter.sh <output_dir> [n_top=30]
#   ROD_FASTA=/path/to/rod.fasta bash diagnose_amf_filter.sh ...   # also compare against ROD
#
# Run in an interactive session, e.g.: srun -N 1 --ntasks=8 --mem=32gb --time=1:00:00 -p interactive --pty bash
# Writes only to <output_dir>/05_asv2otu/04b_amf_filter/diagnostics/

if [[ $# -lt 1 ]]; then
    echo "Usage: bash diagnose_amf_filter.sh <output_dir> [n_top=30]"
    exit 1
fi

OUTPUT_DIR="$(cd "$1" && pwd)"
N_TOP="${2:-30}"
THREADS="${SLURM_NTASKS:-4}"
TARGET_REGEX="${TARGET_REGEX:-Glomeromycot}"

# Keep paths in sync with the Snakefile (SILVA_BLAST_DB, TAXONOMY_DB_PATHS)
TAX_DIR="/projects/standard/kennedyp/shared/taxonomy"
SILVA_DB="${SILVA_DB:-$TAX_DIR/silva_138.2_NR99_blastdb/SILVA_138.2_NR99}"
MAARJAM_FASTA="${MAARJAM_FASTA:-$TAX_DIR/maarjam_dada2_SSU_2021.fasta}"
EUKARYOME_FASTA="${EUKARYOME_FASTA:-$TAX_DIR/DADA2_EUK_SSU_v2.0.fasta}"
ROD_FASTA="${ROD_FASTA:-}"

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"

AMF_DIR="$OUTPUT_DIR/05_asv2otu/04b_amf_filter"
DECISIONS="$AMF_DIR/amf_filter_otu_decisions.tsv"
CURATED_FASTA="$OUTPUT_DIR/05_asv2otu/04_mumu/Centroid_mumu_curated.fas"
DIAG_DIR="$AMF_DIR/diagnostics"

for f in "$DECISIONS" "$CURATED_FASTA"; do
    [[ -f "$f" ]] || { echo "ERROR: not found: $f (has the AMF filter step run?)"; exit 1; }
done

set +u
module load blast-plus/2.13.0-gcc-8.2.0-vo4mr4d 2>/dev/null || true
if ! command -v Rscript &>/dev/null; then
    module load miniforge 2>/dev/null || true
    source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate dada2_processing_v2
fi
set -u
command -v blastn &>/dev/null || { echo "ERROR: blastn not found (module load blast-plus)"; exit 1; }
command -v Rscript &>/dev/null || { echo "ERROR: Rscript not found (conda activate dada2_processing_v2)"; exit 1; }

mkdir -p "$DIAG_DIR"
rm -f "$DIAG_DIR"/blast_*.tsv
TMP_DB_DIR="$(mktemp -d "${TMPDIR:-/tmp}/amf_diag.XXXXXX")"
trap 'rm -rf "$TMP_DB_DIR"' EXIT

echo "AMF filter diagnostics"
echo "  Output dir: $OUTPUT_DIR"
echo "  Top OTUs:   $N_TOP per discard reason"
echo "  Threads:    $THREADS"

# 1. Per-reason summary and query fastas of the most abundant discarded OTUs
Rscript "$PIPELINE_DIR/workflow/scripts/amf_filter_diagnostics.R" select \
    "$DECISIONS" "$CURATED_FASTA" "$DIAG_DIR" "$N_TOP"

# 2. Reference databases: SILVA (prebuilt), others built from fasta into a temp dir
declare -A DBS=()
if [[ -f "${SILVA_DB}.nsq" || -f "${SILVA_DB}.nal" ]]; then
    DBS[SILVA]="$SILVA_DB"
else
    echo "WARNING: SILVA BLAST db not found at $SILVA_DB; skipping"
fi

build_db() {  # <name> <fasta>
    local name="$1" fasta="$2"
    if [[ -z "$fasta" ]]; then
        return
    fi
    if [[ ! -f "$fasta" ]]; then
        echo "WARNING: $name fasta not found at $fasta; skipping"
        return
    fi
    echo "Building temporary $name BLAST db..."
    # U -> T on sequence lines, in case the reference is RNA
    awk '/^>/ {print; next} {gsub(/U/, "T"); gsub(/u/, "t"); print}' "$fasta" > "$TMP_DB_DIR/$name.fasta"
    makeblastdb -in "$TMP_DB_DIR/$name.fasta" -dbtype nucl -out "$TMP_DB_DIR/$name" > "$TMP_DB_DIR/$name.makeblastdb.log"
    rm -f "$TMP_DB_DIR/$name.fasta"
    DBS[$name]="$TMP_DB_DIR/$name"
}
build_db MaarjAM "$MAARJAM_FASTA"
build_db Eukaryome "$EUKARYOME_FASTA"
build_db ROD "$ROD_FASTA"

if [[ ${#DBS[@]} -eq 0 ]]; then
    echo "ERROR: no reference databases available"
    exit 1
fi

# 3. Unfiltered BLAST (no identity/coverage cutoffs; dc-megablast tolerates divergence)
for reason in no_confident_hit non_AMF; do
    query="$DIAG_DIR/top_${reason}.fasta"
    [[ -s "$query" ]] || continue
    for db in "${!DBS[@]}"; do
        echo "BLAST: top $reason OTUs vs $db..."
        blastn -task dc-megablast \
            -query "$query" \
            -db "${DBS[$db]}" \
            -outfmt '6 qseqid sseqid pident length qlen qcovhsp bitscore evalue stitle' \
            -max_target_seqs 20 \
            -max_hsps 1 \
            -evalue 1e-5 \
            -num_threads "$THREADS" \
            -out "$DIAG_DIR/blast_${reason}_${db}.tsv"
    done
done

# 4. Best hit per OTU per database
Rscript "$PIPELINE_DIR/workflow/scripts/amf_filter_diagnostics.R" report \
    "$DECISIONS" "$DIAG_DIR" "$TARGET_REGEX" | tee "$DIAG_DIR/diagnostics_report.txt"

echo ""
echo "Done. Diagnostics in: $DIAG_DIR"
echo "  reason_summary.tsv, length_bins_by_reason.tsv, top_*_besthits.tsv, diagnostics_report.txt"
