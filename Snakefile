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

# Primer sequences (forward and reverse) per amplicon type
PRIMERS = {
    "16S-V4": {
        "fwd": "GTGCCAGCMGCCGCGGTAA",
        "rev": "GGACTACHVGGGTWTCTAAT",
    },
    "ITS1": {
        "fwd": "CTTGGTCATTTAGAGGAAGTAA",
        "rev": "GCTGCGTTCTTCATCGATGC",
    },
    "ITS2": {
        "fwd": "AACTTTYRRCAAYGGATCWCT",
        "rev": "AGCCTCCGCTTATTGATATGCTTAART",
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

# ============================================================================
# Rules
# ============================================================================

rule all:
    """Final target: DADA2 processing complete"""
    input:
        OUTPUT_DIR + f"/03_dada2/{AMPLICON}_combined_sequences_taxa.txt",
        OUTPUT_DIR + "/03_dada2/README.txt"

include: "workflow/rules/cutadapt.smk"
include: "workflow/rules/dada2.smk"
