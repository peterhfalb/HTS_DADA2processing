"""
DADA2 rule: denoising, merging, chimera removal, taxonomy assignment
"""

# ============================================================================
# Helper functions for platform-specific inputs
# ============================================================================

def get_adapter_logs(wildcards):
    """Return adapter trimming logs based on platform"""
    if PLATFORM == "aviti":
        # Aviti uses Trimmomatic logs
        return expand(OUTPUT_DIR + "/01_adapter/01_logs/{sample}_trimmomatic_pass1.log.txt", sample=SAMPLES)
    else:
        # Illumina uses cutadapt logs
        return expand(OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.log.txt", sample=SAMPLES)

def get_adapter_jsons(wildcards):
    """Return adapter trimming JSON logs (Illumina only)"""
    if PLATFORM == "aviti":
        # Aviti doesn't have JSON logs, return empty list (handled by rule)
        return []
    else:
        return expand(OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.json", sample=SAMPLES)

def get_primer_logs(wildcards):
    """Return primer trimming logs based on platform"""
    if PLATFORM == "aviti":
        # For Aviti, primer logs are from cutadapt (middle step)
        return expand(OUTPUT_DIR + "/01b_primer_trimmed/02_logs/cutadapt.{sample}.log.txt", sample=SAMPLES)
    else:
        # For Illumina, primer logs from final trimming
        return expand(OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.log.txt", sample=SAMPLES)

def get_primer_jsons(wildcards):
    """Return primer trimming JSON logs based on platform"""
    if PLATFORM == "aviti":
        # For Aviti, primer JSONs from cutadapt (middle step)
        return expand(OUTPUT_DIR + "/01b_primer_trimmed/02_logs/cutadapt.{sample}.json", sample=SAMPLES)
    else:
        # For Illumina, primer JSONs from final trimming
        return expand(OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.json", sample=SAMPLES)

rule dada2:
    """DADA2 denoising and taxonomy assignment"""
    input:
        reads_r1=expand(
            OUTPUT_DIR + "/02_primer_trimmed/{sample}_R1_001.fastq.gz",
            sample=SAMPLES,
        ),
        reads_r2=expand(
            OUTPUT_DIR + "/02_primer_trimmed/{sample}_R2_001.fastq.gz",
            sample=SAMPLES,
        ),
    output:
        seqtab=OUTPUT_DIR + "/03_dada2/seqtab.rds",
        seqtab_nochim=OUTPUT_DIR + "/03_dada2/seqtab_nochim.rds",
        fasta=OUTPUT_DIR + "/03_dada2/sequences.fasta",
        summary=OUTPUT_DIR + "/03_dada2/sequence_process_summary.txt",
        taxa=OUTPUT_DIR + "/03_dada2/taxID.rds",
        taxa_boot=OUTPUT_DIR + "/03_dada2/taxID_bootstrap.rds",
        errF=OUTPUT_DIR + "/03_dada2/errF.rds",
        errR=OUTPUT_DIR + "/03_dada2/errR.rds",
        combined_taxa=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        combined_taxa_boot=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_bootstrap_{DB_NAME}.txt",
    params:
        outdir=OUTPUT_DIR + "/03_dada2",
        amplicon=AMPLICON,
        quality=QUALITY,
        taxonomy_db=EFFECTIVE_TAXONOMY_DB,
        project_name=PROJECT_NAME,
        db_name=DB_NAME,
        platform=PLATFORM,
        fwd_reads_only=int(FWD_READS_ONLY),
    log:
        OUTPUT_DIR + "/.logs/dada2.log",
    threads: workflow.cores
    shell:
        """
        mkdir -p {params.outdir}
        Rscript workflow/scripts/run_dada2.R \
          {params.outdir} \
          {params.amplicon} \
          {params.quality} \
          {params.taxonomy_db} \
          {threads} \
          {params.project_name} \
          {params.db_name} \
          {params.platform} \
          {params.fwd_reads_only} \
          2>&1 | tee {log}
        """

rule dada2_readme:
    """Create README for DADA2 output directory"""
    output:
        OUTPUT_DIR + "/03_dada2/README.txt",
    params:
        amplicon=AMPLICON,
        quality=QUALITY,
        project_name=PROJECT_NAME,
        db_name=DB_NAME,
    shell:
        """
        cat > {output} << EOF
# DADA2 Denoising and Taxonomy Assignment Results
## Directory: 03_dada2/

This directory contains the final DADA2 output including denoised ASVs and taxonomy assignments.

### Contents:

**Sequence Tables (R objects):**
- `seqtab.rds`: Sequence table before chimera removal
- `seqtab_nochim.rds`: Sequence table after chimera removal (main output)

**Fasta:**
- `sequences.fasta`: Unique sequences from seqtab_nochim in fasta format

**Summary Statistics:**
- `sequence_process_summary.txt`: Per-sample counts at each processing step
  Columns: dada2_input, filtered, dada_f, dada_r, merged, nonchim

**Taxonomy (R objects):**
- `taxID.rds`: Taxonomy assignments
- `taxID_bootstrap.rds`: Bootstrap confidence values for taxonomy

**Combined Results:**
- `{params.project_name}__combined_sequences_ASVtaxa_{params.db_name}.txt`: ASV abundance table with taxonomy
- `{params.project_name}__combined_sequences_ASVtaxa_bootstrap_{params.db_name}.txt`: ASV abundance with taxonomy and bootstrap values

**Processing Details:**
- Amplicon type: {params.amplicon}
- Quality setting: {params.quality}
- Filtering: Based on {params.quality} quality parameters
- Dereplication: Yes
- Denoising: DADA2 with pseudo-pooling
- Merging: Paired reads merged with minimum overlap and trimming
- Chimera removal: Bimera, consensus method
- Taxonomy database: MSI Agate shared location (pr2/silva/unite/maarjam)

### Next Steps:
Use the combined_sequences_taxa.txt files for downstream analysis (diversity, relative abundance, etc.)
EOF
        """

rule dada2_qc:
    """Generate QC stats file and figures"""
    input:
        dada2_summary = OUTPUT_DIR + "/03_dada2/sequence_process_summary.txt",
        errF          = OUTPUT_DIR + "/03_dada2/errF.rds",
        errR          = OUTPUT_DIR + "/03_dada2/errR.rds",
        seqtab_nochim = OUTPUT_DIR + "/03_dada2/seqtab_nochim.rds",
        adapter_jsons = get_adapter_jsons,
        primer_jsons  = get_primer_jsons,
        adapter_logs  = get_adapter_logs,
        primer_logs   = get_primer_logs,
        primer_fqs    = expand(OUTPUT_DIR + "/02_primer_trimmed/{sample}_R1_001.fastq.gz", sample=SAMPLES),
    output:
        qc_summary = OUTPUT_DIR + "/04_dada2_QCsummary/qc_summary.txt",
    params:
        outdir    = OUTPUT_DIR + "/04_dada2_QCsummary",
        amplicon  = AMPLICON,
        platform  = PLATFORM,
        dada2_dir = OUTPUT_DIR + "/03_dada2",
        primer_trimmed_dir = OUTPUT_DIR + "/02_primer_trimmed",
    log:
        OUTPUT_DIR + "/.logs/qc.log",
    shell:
        """
        mkdir -p {params.outdir}/figures
        Rscript workflow/scripts/qc_figures.R \
          {params.outdir} {params.amplicon} {params.platform} \
          {input.dada2_summary} {params.dada2_dir} \
          {params.primer_trimmed_dir} \
          2>&1 | tee {log}
        """
