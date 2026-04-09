#!/bin/bash
set -euo pipefail

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

# Summary
echo ""
echo "=========================================="
echo "✓ SETUP COMPLETE"
echo "=========================================="
echo ""
echo "The environment is ready. You can now submit SLURM jobs:"
echo ""
echo "  cd $PIPELINE_DIR"
echo "  run_dada2processing --amplicon 16S-V4 --fastq-dir ./raw_fastqs \\"
echo "    --output-dir ./output --email user@example.com"
echo ""
echo "Then:"
echo "  sbatch ./output/submit_job.sh"
echo ""
