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
OUTPUT_DIR = config["output_dir"]
SAMPLES = config["samples"]
EMAIL = config["email"]
FASTQ_DIR = config["fastq_dir"]
PROJECT_NAME = config.get("project_name", os.path.basename(OUTPUT_DIR))
TAXONOMY_DB_OVERRIDE = config.get("taxonomy_db_override", "")
FWD_PRIMER_OVERRIDE = config.get("fwd_primer_override", "")
REV_PRIMER_OVERRIDE = config.get("rev_primer_override", "")

# Primer sequences (forward and reverse) per amplicon type
PRIMERS = {
    "16S-V4": {
        "fwd": "GTGYCAGCMGCCGCGGTAA",
        "rev": "GGACTACNVGGGTWTCTAAT",
    },
    "ITS1": {
        "fwd": "CTTGGTCATTTAGAGGAAGTAA",
        "rev": "GCTGCGTTCTTCATCGATGC",
    },
    "ITS2": {
        "fwd": "TCGATGAAGAACGCAGCG",
        "rev": "TCCTCCGCTTATTGATATGC",
    },
    "18S-AMF": {
        "fwd": "CAGCCGCGGTAATTCCAGCT",
        "rev": "GAACCCAAACACTTTGGTTTCC",
    },
    "18S-V4": {
        "fwd": "TTAAARVGYTCGTAGTYG",
        "rev": "CCGTCAATTHCTTYAART",
    },
}

# Taxonomy database paths on MSI Agate
TAXONOMY_DB = {
    "16S-V4": "/projects/standard/kennedyp/shared/taxonomy/silva_nr99_v138.1_train_set.fa",
    "ITS1": "/projects/standard/kennedyp/shared/taxonomy/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta",
    "ITS2": "/projects/standard/kennedyp/shared/taxonomy/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta",
    "18S-AMF": "/projects/standard/kennedyp/shared/taxonomy/maarjam_dada2_SSU_2021.fasta",
    "18S-V4": "/projects/standard/kennedyp/shared/taxonomy/pr2_version_5.1.1_SSU_dada2.fasta",
}

# Named taxonomy database paths (for --taxonomy-database option)
TAXONOMY_DB_PATHS = {
    "SILVA": "/projects/standard/kennedyp/shared/taxonomy/silva_nr99_v138.1_train_set.fa",
    "UNITE": "/projects/standard/kennedyp/shared/taxonomy/sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta",
    "PR2": "/projects/standard/kennedyp/shared/taxonomy/pr2_version_5.1.1_SSU_dada2.fasta",
    "MaarjAM": "/projects/standard/kennedyp/shared/taxonomy/maarjam_dada2_SSU_2021.fasta",
    "EukaryomeSSU": "/projects/standard/kennedyp/shared/taxonomy/DADA2_EUK_SSU_v.2.0.fasta",
    "EukaryomeITS": "/projects/standard/kennedyp/shared/taxonomy/DADA2_EUK_ITS_v.2.0_wKennedySynmock.fasta",
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
EFFECTIVE_TAXONOMY_DB = TAXONOMY_DB_PATHS[TAXONOMY_DB_OVERRIDE] if TAXONOMY_DB_OVERRIDE else TAXONOMY_DB[AMPLICON]
DB_NAME = TAXONOMY_DB_OVERRIDE if TAXONOMY_DB_OVERRIDE else TAXONOMY_DB_DEFAULT_NAMES[AMPLICON]

# Resolve effective primers
EFFECTIVE_FWD = FWD_PRIMER_OVERRIDE if FWD_PRIMER_OVERRIDE else PRIMERS[AMPLICON]["fwd"]
EFFECTIVE_REV = REV_PRIMER_OVERRIDE if REV_PRIMER_OVERRIDE else PRIMERS[AMPLICON]["rev"]

# ============================================================================
# Rules
# ============================================================================

rule all:
    """Final target: DADA2 processing complete"""
    input:
        OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_bootstrap_{DB_NAME}.txt",
        OUTPUT_DIR + "/03_dada2/README.txt",
        OUTPUT_DIR + "/04_QC/qc_summary.txt"

include: "workflow/rules/cutadapt.smk"
include: "workflow/rules/dada2.smk"
