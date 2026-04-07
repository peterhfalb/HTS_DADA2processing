"""
DADA2 rule: denoising, merging, chimera removal, taxonomy assignment
"""

rule dada2:
    """DADA2 denoising and taxonomy assignment"""
    input:
        reads=expand(
            OUTPUT_DIR + "/02_primer_trimmed/{sample}_L001_R1_001.fastq.gz",
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
          2>&1 | tee {log}
        """

rule dada2_readme:
    """Create README for DADA2 output directory"""
    output:
        OUTPUT_DIR + "/03_dada2/README.txt",
    params:
        amplicon=AMPLICON,
        quality=QUALITY,
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
- `{params.amplicon}_combined_sequences_taxa.txt`: Sequence abundance table with taxonomy
- `{params.amplicon}_combined_sequences_taxa_bootstrap.txt`: Sequence abundance with taxonomy and bootstrap values

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
        adapter_logs  = expand(OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.log.txt", sample=SAMPLES),
        primer_logs   = expand(OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.log.txt", sample=SAMPLES),
        primer_fqs    = expand(OUTPUT_DIR + "/02_primer_trimmed/{sample}_L001_R1_001.fastq.gz", sample=SAMPLES),
    output:
        qc_summary = OUTPUT_DIR + "/04_QC/qc_summary.txt",
    params:
        outdir    = OUTPUT_DIR + "/04_QC",
        amplicon  = AMPLICON,
        dada2_dir = OUTPUT_DIR + "/03_dada2",
        primer_trimmed_dir = OUTPUT_DIR + "/02_primer_trimmed",
    log:
        OUTPUT_DIR + "/.logs/qc.log",
    shell:
        """
        mkdir -p {params.outdir}/figures
        Rscript workflow/scripts/qc_figures.R \
          {params.outdir} {params.amplicon} \
          {input.dada2_summary} {params.dada2_dir} \
          {params.primer_trimmed_dir} \
          2>&1 | tee {log}
        """
