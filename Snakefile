"""
DADA2 Snakemake Pipeline
Amplicon denoising and taxonomy assignment for 16S-V4, ITS1, ITS2, 18S-AMF, 18S-V4
"""

import os
from pathlib import Path
from datetime import datetime

# ============================================================================
# Configuration
# ============================================================================

configfile: "config.yaml"

AMPLICON = config["amplicon"]
QUALITY = config["quality"]
PLATFORM = config["platform"]
OUTPUT_DIR = config["output_dir"]
SAMPLES = config["samples"]
EMAIL = config["email"]
FASTQ_DIR = config["fastq_dir"]
PROJECT_NAME = config["project_name"]
TAXONOMY_DB_OVERRIDE = config.get("taxonomy_db_override", "")
FWD_PRIMER_OVERRIDE = config.get("fwd_primer_override", "")
REV_PRIMER_OVERRIDE = config.get("rev_primer_override", "")

# DADA2 processing configuration (handle both string and int YAML inputs)
FWD_READS_ONLY  = str(config.get("fwd_reads_only", "0")) == "1"

# OTU pipeline configuration (handle both string and int YAML inputs)
SKIP_OTU        = str(config.get("skip_otu", "0")) == "1"
RUN_ITSX        = str(config.get("run_itsx", "0")) == "1"
OTU_CLUSTER_ID  = config.get("otu_cluster_id", "0.97")
MUMU_BLAST_ID   = config.get("mumu_blast_id", "84")
MUMU_RATIO      = config.get("mumu_ratio", "1")

# Load primer sequences from external config file (single source of truth)
import json
primers_config_path = os.path.join(workflow.basedir, "workflow/config/primers.json")
with open(primers_config_path) as f:
    PRIMERS = json.load(f)

# Named taxonomy database paths (for --taxonomy-database option)
TAXONOMY_DB_PATHS = {
    "SILVA": "/projects/standard/kennedyp/shared/taxonomy/silva_nr99_v138.1_train_set.fa",
    "UNITE": "/projects/standard/kennedyp/shared/taxonomy/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta",
    "PR2": "/projects/standard/kennedyp/shared/taxonomy/pr2_version_5.1.1_SSU_dada2.fasta",
    "MaarjAM": "/projects/standard/kennedyp/shared/taxonomy/maarjam_dada2_SSU_2021.fasta",
    "EukaryomeSSU": "/projects/standard/kennedyp/shared/taxonomy/DADA2_EUK_SSU_v2.0.fasta",
    "EukaryomeITS": "/projects/standard/kennedyp/shared/taxonomy/DADA2_EUK_ITS_v2.0_wKennedySynmock.fasta",
}

# Default DB name per amplicon (for output filename when no override)
TAXONOMY_DB_DEFAULT_NAMES = {
    "16S-V4": "SILVA",
    "ITS1": "UNITE",
    "ITS2": "UNITE",
    "18S-AMF": "MaarjAM",
    "18S-V4": "PR2",
}

# Resolve effective taxonomy database path and short name
if TAXONOMY_DB_OVERRIDE:
    if TAXONOMY_DB_OVERRIDE not in TAXONOMY_DB_PATHS:
        raise KeyError(f"taxonomy_db_override '{TAXONOMY_DB_OVERRIDE}' not in TAXONOMY_DB_PATHS. Valid keys: {list(TAXONOMY_DB_PATHS.keys())}")
    EFFECTIVE_TAXONOMY_DB = TAXONOMY_DB_PATHS[TAXONOMY_DB_OVERRIDE]
    DB_NAME = TAXONOMY_DB_OVERRIDE
else:
    if AMPLICON not in TAXONOMY_DB_DEFAULT_NAMES:
        raise KeyError(f"amplicon '{AMPLICON}' not in TAXONOMY_DB_DEFAULT_NAMES. Valid amplicons: {list(TAXONOMY_DB_DEFAULT_NAMES.keys())}")
    db_name_key = TAXONOMY_DB_DEFAULT_NAMES[AMPLICON]
    EFFECTIVE_TAXONOMY_DB = TAXONOMY_DB_PATHS[db_name_key]
    DB_NAME = db_name_key

# Resolve effective primers
EFFECTIVE_FWD = FWD_PRIMER_OVERRIDE if FWD_PRIMER_OVERRIDE else PRIMERS[AMPLICON]["fwd"]
EFFECTIVE_REV = REV_PRIMER_OVERRIDE if REV_PRIMER_OVERRIDE else PRIMERS[AMPLICON]["rev"]

# Load detected adapters from adapter auto-detection
# Note: includes transposase sequences for read-through contamination on short amplicons (16S-V4, ITS)
DETECTED_ADAPTERS = {
    "adapters": [
        "CTGTCTCTTATACACATCT",           # Nextera transposase (read-through R1, Illumina)
        "AGATGTGTATAAGAGACAG",           # Nextera transposase RC (read-through R2, Illumina)
        "ATCTCGTATGCCGTCTTCTGCTTG",      # i7 Nextera XT index adapter
        "CAAGCAGAAGACGGCATACGAGAT",      # i7 index adapter RC
        "GTGTAGATCTCGGTGGTCGCCGTATCATT", # i5 Nextera XT index adapter
        "AATGATACGGCGACCACCGAGATCTACAC", # i5 index adapter RC
    ]
}

# Try to load auto-detected adapters if available
adapter_detection_json = os.path.join(OUTPUT_DIR, ".adapter_detection.json")
if os.path.exists(adapter_detection_json):
    try:
        import json
        with open(adapter_detection_json) as f:
            adapter_data = json.load(f)
            # Use the final_adapters list which combines detected + defaults
            if "final_adapters" in adapter_data and adapter_data["final_adapters"]:
                DETECTED_ADAPTERS["adapters"] = adapter_data["final_adapters"]
    except Exception as e:
        print(f"Warning: Could not load adapter detection results: {e}")
        print("Using default Nextera XT adapters")

# ============================================================================
# Helper Functions
# ============================================================================

def get_input_r1(wildcards):
    """Find R1 file with or without _L001_ lane number"""
    sample = wildcards.sample
    with_lane = OUTPUT_DIR + f"/00_raw/{sample}_L001_R1_001.fastq.gz"
    without_lane = OUTPUT_DIR + f"/00_raw/{sample}_R1_001.fastq.gz"
    if os.path.exists(with_lane):
        return with_lane
    elif os.path.exists(without_lane):
        return without_lane
    else:
        return with_lane  # Return default; let Snakemake error

def get_input_r2(wildcards):
    """Find R2 file with or without _L001_ lane number"""
    sample = wildcards.sample
    with_lane = OUTPUT_DIR + f"/00_raw/{sample}_L001_R2_001.fastq.gz"
    without_lane = OUTPUT_DIR + f"/00_raw/{sample}_R2_001.fastq.gz"
    if os.path.exists(with_lane):
        return with_lane
    elif os.path.exists(without_lane):
        return without_lane
    else:
        return with_lane  # Return default; let Snakemake error

# ============================================================================
# Rules
# ============================================================================

rule cleanup_temp_files:
    """Remove temporary FASTQ files created during processing"""
    input:
        # ASV outputs (always)
        OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_bootstrap_{DB_NAME}.txt",
        OUTPUT_DIR + "/03_dada2/README.txt",
        OUTPUT_DIR + "/04_dada2_QCsummary/qc_summary.txt",
        # Main directory summary files (always)
        OUTPUT_DIR + "/qc_summary.txt",
        OUTPUT_DIR + "/README.txt"
    output:
        touch_file=OUTPUT_DIR + "/.cleanup.done"
    shell:
        """
        echo "Cleaning up temporary files..."
        rm -rf {OUTPUT_DIR}/01_adapter
        rm -rf {OUTPUT_DIR}/01b_primer_trimmed
        rm -rf {OUTPUT_DIR}/02_primer_trimmed
        rm -rf {OUTPUT_DIR}/03_dada2/filtered
        echo "Temporary file cleanup complete"
        touch {output.touch_file}
        """

rule all:
    """Final target: DADA2 + OTU processing complete"""
    input:
        OUTPUT_DIR + "/.cleanup.done"

# Include trimming rules based on platform
if PLATFORM == "aviti":
    include: "workflow/rules/aviti_trimming.smk"
else:
    include: "workflow/rules/cutadapt.smk"

include: "workflow/rules/dada2.smk"
include: "workflow/rules/asv2otu.smk"
