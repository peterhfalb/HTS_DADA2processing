#!/bin/bash
set -euo pipefail

RED='\033[0;31m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
WHITE='\033[1;37m'
RESET='\033[0m'

WIDTH=28
STEPS=60

echo ""

for ((i=0; i<STEPS; i++)); do
    left=$(awk -v i="$i" -v w="$WIDTH" 'BEGIN { printf "%d", int((sin(i * 0.27) + 1) / 2 * w) }')
    right=$(awk -v i="$i" -v w="$WIDTH" 'BEGIN { printf "%d", int((sin(i * 0.27 + 3.14159) + 1) / 2 * w) }')

    if [ "$left" -gt "$right" ]; then
        tmp=$left; left=$right; right=$tmp
    fi

    line=""
    for ((col=0; col<=WIDTH; col++)); do
        if [ "$col" -eq "$left" ] && [ "$col" -eq "$right" ]; then
            line+="${YELLOW}X${RESET}"
        elif [ "$col" -eq "$left" ]; then
            line+="${RED}O${RESET}"
        elif [ "$col" -eq "$right" ]; then
            line+="${CYAN}O${RESET}"
        elif [ "$col" -gt "$left" ] && [ "$col" -lt "$right" ]; then
            if (( col % 2 == 0 )); then
                line+="${YELLOW}-${RESET}"
            else
                line+="${WHITE}-${RESET}"
            fi
        else
            line+=" "
        fi
    done

    echo -e "                  $line"
    sleep 0.03
done

echo ""

# ── Banner ─────────────────────────────────────────────────────────────────────
# Visual box width = 64: '  ║' + 60 chars inner + '║'
# Content rows: '  ║  ' (5) + art + padding + '║'  → art+padding = 58
# DADA2 art = 40 wide  → 18 spaces padding
# PIPELINE art = 55 wide → 3 spaces padding
# Tagline = 54 wide (🧬 counts as 2 each) → 3 left + 3 right padding
BRED='\033[1;31m'
BYELLOW='\033[1;33m'
BCYAN='\033[1;36m'

bl() { echo -e "$1"; sleep 0.04; }

sleep 0.15

E='                  '  # 18 spaces (DADA2 padding)
P='   '                 # 3 spaces  (PIPELINE padding)

bl "${BCYAN}  ╔════════════════════════════════════════════════════════════╗${RESET}"
bl "${BCYAN}  ║${RESET}                                                            ${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}██████╗  █████╗ ██████╗  █████╗ ██████╗ ${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}██╔══██╗██╔══██╗██╔══██╗██╔══██╗╚════██╗${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}██║  ██║███████║██║  ██║███████║ █████╔╝${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}██║  ██║██╔══██║██║  ██║██╔══██║██╔═══╝ ${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}██████╔╝██║  ██║██████╔╝██║  ██║███████╗${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BRED}╚═════╝ ╚═╝  ╚═╝╚═════╝ ╚═╝  ╚═╝╚══════╝${RESET}${E}${BCYAN}║${RESET}"
bl "${BCYAN}  ║${RESET}                                                            ${BCYAN}║${RESET}"
bl "${BCYAN}  ╠════════════════════════════════════════════════════════════╣${RESET}"
bl "${BCYAN}  ║${RESET}                                                            ${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}██████╗ ██╗██████╗ ███████╗██╗     ██╗███╗  ██╗███████╗${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}██╔══██╗██║██╔══██╗██╔════╝██║     ██║████╗ ██║██╔════╝${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}██████╔╝██║██████╔╝█████╗  ██║     ██║██╔██╗██║█████╗  ${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}██╔═══╝ ██║██╔═══╝ ██╔══╝  ██║     ██║██║╚████║██╔══╝  ${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}██║     ██║██║     ███████╗███████╗██║██║ ╚███║███████╗${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║  ${BYELLOW}╚═╝     ╚═╝╚═╝     ╚══════╝╚══════╝╚═╝╚═╝  ╚══╝╚══════╝${RESET}${P}${BCYAN}║${RESET}"
bl "${BCYAN}  ║${RESET}                                                            ${BCYAN}║${RESET}"
bl "${BCYAN}  ╠════════════════════════════════════════════════════════════╣${RESET}"
sleep 0.1
bl "${BCYAN}  ║${RESET}   🧬  Amplicon denoising - It's a SLURM-apocalypse!!  🧬   ${BCYAN}║${RESET}"
bl "${BCYAN}  ╚════════════════════════════════════════════════════════════╝${RESET}"
echo ""


# DADA2 Pipeline Installation Script
# Sets up conda environment with all required dependencies

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_FILE="$PIPELINE_DIR/envs/dada2.yaml"
ENV_NAME="dada2_processing_v2"

echo "=========================================="
echo "DADA2 Pipeline Environment Setup"
echo "=========================================="
echo ""

# Load miniforge (modern conda/mamba)
module load miniforge/24.3

# Source conda/mamba
CONDA_BASE=$(conda info --base 2>/dev/null) || { echo "ERROR: miniforge not found."; exit 1; }
source "$CONDA_BASE/etc/profile.d/conda.sh"

# Check if environment file exists
if [[ ! -f "$ENV_FILE" ]]; then
    echo "ERROR: Environment file not found: $ENV_FILE"
    exit 1
fi

echo "Pipeline directory: $PIPELINE_DIR"
echo "Environment file:   $ENV_FILE"
echo "Environment name:   $ENV_NAME"
echo ""

# Check if environment already exists
if conda env list | grep -q "^$ENV_NAME "; then
    echo "⚠️  Environment '$ENV_NAME' already exists"
    read -p "Delete and recreate? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo "Removing existing environment..."
        conda env remove -n "$ENV_NAME" --yes
    else
        echo "Using existing environment. Skipping creation."
        SKIP_CREATE=1
    fi
fi

# Create environment if needed
if [[ "${SKIP_CREATE:-0}" != "1" ]]; then
    echo "Creating conda environment (this may take 5-15 minutes)..."
    echo ""
    conda env create -f "$ENV_FILE" --yes
    echo ""
    echo "✓ Environment created"
fi

# ==============================================================================
# VERIFICATION
# ==============================================================================
echo ""
echo "=========================================="
echo "VERIFYING PACKAGES"
echo "=========================================="
echo ""

# Activate environment
echo "Activating environment..."
conda activate "$ENV_NAME"

# Verify system modules (outside conda environment)
echo ""
echo "Checking system modules..."

set +u
module load cutadapt/5.0 2>/dev/null && {
    CUTADAPT_VERSION=$(cutadapt --version 2>&1 | head -1)
    echo "✓ cutadapt/5.0 ($CUTADAPT_VERSION)"
} || {
    echo "⚠️  cutadapt/5.0 module not available (will load in SLURM jobs)"
}

module load snakemake 2>/dev/null && {
    SNAKEMAKE_VERSION=$(snakemake --version 2>&1)
    echo "✓ snakemake ($SNAKEMAKE_VERSION)"
} || {
    echo "⚠️  snakemake module not available (will load in SLURM jobs)"
}
set -u

# Verify R packages
echo ""
echo "Checking R packages..."

R_VERIFY=$(R --quiet --slave << 'RSCRIPT'
packages <- c("dada2", "jsonlite", "ggplot2")
results <- list()
for (pkg in packages) {
    results[[pkg]] <- require(pkg, character.only = TRUE, quietly = TRUE)
}
for (pkg in packages) {
    status <- ifelse(results[[pkg]], "✓", "✗")
    cat(sprintf("%s %s\n", status, pkg))
    if (!results[[pkg]]) {
        quit(status = 1)
    }
}
RSCRIPT
)

if [[ $? -eq 0 ]]; then
    echo "$R_VERIFY"
else
    echo "✗ Some R packages failed to load"
    echo "$R_VERIFY"
    exit 1
fi

# ==============================================================================
# SETTING UP MUMU
# ==============================================================================
echo ""
echo "=========================================="
echo "SETTING UP MUMU"
echo "=========================================="
echo ""

MUMU_DIR="$HOME/packages/mumu"
MUMU_BIN="$MUMU_DIR/mumu"

if [[ -f "$MUMU_BIN" ]]; then
    echo "✓ mumu already compiled at $MUMU_BIN"
else
    echo "Compiling mumu from source..."
    set +u
    module load gcc/13.1.0-5z64cho 2>/dev/null || echo "WARNING: gcc/13.1.0 not available"
    set -u
    mkdir -p "$(dirname "$MUMU_DIR")"
    git clone https://github.com/frederic-mahe/mumu.git "$MUMU_DIR" 2>&1
    cd "$MUMU_DIR" && make 2>&1
    cd "$PIPELINE_DIR"
    if [[ -f "$MUMU_BIN" ]]; then
        echo "✓ mumu compiled at $MUMU_BIN"
    else
        echo "✗ mumu compilation failed"
        echo "  You can try manually with:"
        echo "  cd $MUMU_DIR && make"
        exit 1
    fi
fi

# Summary
echo ""
echo "=========================================="
echo "✓ SETUP COMPLETE"
echo "=========================================="
echo ""

# ==============================================================================
# ADD COMMAND TO PATH
# ==============================================================================
BIN_DIR="$PIPELINE_DIR/bin"
BASHRC="$HOME/.bashrc"

if [[ ! -d "$BIN_DIR" ]]; then
    echo "✗ bin directory not found at $BIN_DIR"
    exit 1
fi

echo "Adding pipeline command to PATH..."
if grep -q "DADA2_PIPELINE_BIN" "$BASHRC"; then
    echo "✓ PATH already configured in bashrc"
else
    echo "" >> "$BASHRC"
    echo "# DADA2 Pipeline" >> "$BASHRC"
    echo "export DADA2_PIPELINE_BIN=\"$BIN_DIR\"" >> "$BASHRC"
    echo "export PATH=\"\$DADA2_PIPELINE_BIN:\$PATH\"" >> "$BASHRC"
    echo "✓ Added pipeline bin directory to PATH in $BASHRC"
fi

echo ""
echo "The environment is ready. You can now submit SLURM jobs!!!"
echo ""
echo "⚠️  To use run_dada2processing, log out and back in (or start a new SSH session)."
echo "    'source ~/.bashrc' may not fully refresh the shell environment."
echo ""
echo ""
echo "Usage: run_dada2processing [OPTIONS]"
echo ""
echo "Required arguments:"
echo "  --amplicon AMPLICON       Amplicon type: 16S-V4, ITS1, ITS2, 18S-AMF, 18S-V4"
echo "  --fastq-dir DIR           Directory containing raw demultiplexed fastq files"
echo "  --output-dir DIR          Output directory (will be created if needed)"
echo "  --email EMAIL             Email for SLURM notifications"
echo "  --platform PLATFORM       Sequencing platform: illumina or aviti"
echo "  --project-name NAME       Project name for output file naming"
echo ""
echo "Optional arguments:"
echo "  --quality QUALITY         Sequence quality level: good (default) or bad"
echo "  --taxonomy-database DB    Override taxonomy database: SILVA, UNITE, PR2, MaarjAM, EukaryomeSSU, EukaryomeITS"
echo "  --fwd-primer SEQUENCE     Override forward primer for this amplicon (DNA string)"
echo "  --rev-primer SEQUENCE     Override reverse primer for this amplicon (DNA string)"
echo "  --fwd-reads-only          Use forward reads only; skip reverse read processing (default: off)"
echo "  --skip-otu                Skip OTU clustering pipeline (default: run OTU pipeline)"
echo "  --run-itsx                Run ITSx extraction (ITS1/ITS2 only; default: off)"
echo "  --otu-cluster-id ID       VSEARCH clustering identity (default: 0.97)"
echo "  --mumu-blast-id PCT       BLAST identity % for mumu (default: 94 for 16S-V4, 84 for others)"
echo "  --mumu-ratio RATIO        mumu minimum_ratio (default: 100 for 16S-V4, 1 for others)"
echo "  --jobs JOBS               Max parallel Snakemake jobs (default: 16)"
echo "  --slurm-time TIME         SLURM walltime (HH:MM:SS, default: 12:00:00, 24:00:00 for 16S-V4 or Aviti)"
echo "  --slurm-mem MEM           SLURM memory allocation (default: 80GB + 20GB per 100 samples, max 240GB)"
echo "  --help                    Show this help message"
echo ""
echo "Examples:"
echo "  run_dada2processing --amplicon 16S-V4 --fastq-dir /mnt/data/raw_fastqs --output-dir /mnt/output/16s_run --email me@umn.edu --platform illumina --project-name MyProject"
echo "  run_dada2processing --amplicon ITS2 --fastq-dir ./fastqs --output-dir ./output --email me@umn.edu --platform aviti --project-name ITS_Study --quality bad"
echo "  run_dada2processing --amplicon 16S-V4 --fastq-dir ./fastqs --output-dir ./output --email me@umn.edu --platform illumina --project-name MyProject --taxonomy-database UNITE"
