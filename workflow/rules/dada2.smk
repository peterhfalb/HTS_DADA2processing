"""
DADA2 rules: denoising, merging, chimera removal, AMF SILVA screen (18S-AMF only), taxonomy assignment
"""

# ============================================================================
# Helper functions for platform-specific inputs
# ============================================================================

def get_adapter_logs(wildcards):
    """Return adapter cutadapt logs (same location on both platforms: trim_adapters on
    Illumina, aviti_adapter_stats on Aviti)"""
    return expand(OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.log.txt", sample=SAMPLES)

def get_adapter_jsons(wildcards):
    """Return adapter cutadapt JSON logs (same location on both platforms)"""
    return expand(OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.json", sample=SAMPLES)

def get_primer_logs(wildcards):
    """Return primer trimming logs based on platform"""
    if PLATFORM == "aviti":
        # For Aviti, primer logs are from cutadapt (middle step), stored in .logs/
        return expand(OUTPUT_DIR + "/.logs/aviti_cutadapt_primers_{sample}.log", sample=SAMPLES)
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

def get_asv_taxonomy_input(wildcards):
    """Sequence table for ASV taxonomy (AMF SILVA screen conditional)"""
    if AMF_FILTER:
        return OUTPUT_DIR + "/03_dada2/amf_filter/seqtab_nochim_AMF_filtered.rds"
    return OUTPUT_DIR + "/03_dada2/seqtab_nochim.rds"

rule dada2:
    """DADA2 denoising"""
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
        errF=OUTPUT_DIR + "/03_dada2/errF.rds",
        errR=OUTPUT_DIR + "/03_dada2/errR.rds",
        asv_counts=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__ASV_counts_unfiltered.txt",
    params:
        outdir=OUTPUT_DIR + "/03_dada2",
        amplicon=AMPLICON,
        quality=QUALITY,
        project_name=PROJECT_NAME,
        platform=PLATFORM,
        fwd_reads_only=int(FWD_READS_ONLY),
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/dada2.log",
    threads: workflow.cores
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.pipeline_dir}/workflow/scripts/run_dada2.R \
          {params.outdir} \
          {params.amplicon} \
          {params.quality} \
          {threads} \
          {params.project_name} \
          {params.platform} \
          {params.fwd_reads_only} \
          2>&1 | tee {log}
        """

# Only define the ASV AMF SILVA screen for 18S-AMF with the filter enabled
if AMF_FILTER:
    rule asv_amf_silva_filter:
        """BLAST ASVs against SILVA; keep ASVs whose near-top hits are mostly AMF (same logic as amf_silva_filter)"""
        input:
            seqtab=OUTPUT_DIR + "/03_dada2/seqtab_nochim.rds",
        output:
            seqtab=OUTPUT_DIR + "/03_dada2/amf_filter/seqtab_nochim_AMF_filtered.rds",
            fasta=OUTPUT_DIR + "/03_dada2/amf_filter/ASVs_AMF_filtered.fasta",
            decisions=OUTPUT_DIR + "/03_dada2/amf_filter/amf_filter_asv_decisions.tsv",
            discarded=OUTPUT_DIR + "/03_dada2/amf_filter/amf_filter_discarded.tsv",
            per_sample=OUTPUT_DIR + "/03_dada2/amf_filter/amf_filter_per_sample.tsv",
            counts=OUTPUT_DIR + "/03_dada2/amf_filter/amf_filter_counts.txt",
        params:
            outdir=OUTPUT_DIR + "/03_dada2/amf_filter",
            project=PROJECT_NAME,
            blast_db=SILVA_BLAST_DB,
            min_id=AMF_MIN_ID,
            min_qcov=AMF_MIN_QCOV,
            top_pct=AMF_TOP_PCT,
            min_frac=AMF_MIN_FRAC,
            target_regex=AMF_TARGET_REGEX,
            pipeline_dir=workflow.basedir,
        threads: workflow.cores
        log:
            OUTPUT_DIR + "/.logs/asv_amf_silva_filter.log",
        shell:
            """
            mkdir -p {params.outdir}

            # Verify SILVA BLAST database exists
            if [[ ! -f "{params.blast_db}.nsq" && ! -f "{params.blast_db}.nal" ]]; then
                echo "ERROR: SILVA BLAST database not found at {params.blast_db}" | tee {log}
                echo "Build it once with: bash {params.pipeline_dir}/bin/setup_silva_blastdb.sh" | tee -a {log}
                exit 1
            fi

            Rscript {params.pipeline_dir}/workflow/scripts/asv_amf_filter.R prepare \
                {input.seqtab} {params.outdir} 2>&1 | tee {log}

            # BLAST against SILVA (hits restricted to identity and query coverage cutoffs)
            blastn -query {params.outdir}/asv_sequences.fasta \
                -db {params.blast_db} \
                -outfmt '6 qseqid sseqid pident length qcovhsp bitscore evalue stitle' \
                -max_target_seqs 100 \
                -max_hsps 1 \
                -perc_identity {params.min_id} \
                -qcov_hsp_perc {params.min_qcov} \
                -num_threads {threads} \
                -out {params.outdir}/silva_blast_hits.tsv 2>&1 | tee -a {log}

            # Decide keep/discard per ASV and write QC records
            Rscript {params.pipeline_dir}/workflow/scripts/amf_silva_filter.R \
                {params.outdir}/silva_blast_hits.tsv \
                {params.outdir}/asv_sequences.fasta \
                {params.outdir}/asv_counts_with_ids.txt \
                {params.outdir}/asv_empty_map.csv \
                {params.outdir} \
                {params.project} \
                {params.top_pct} \
                {params.min_frac} \
                '{params.target_regex}' \
                ASV \
                2>&1 | tee -a {log}

            # Subset the sequence table to kept ASVs
            Rscript {params.pipeline_dir}/workflow/scripts/asv_amf_filter.R finalize \
                {input.seqtab} {params.outdir} {output.seqtab} {params.project} 2>&1 | tee -a {log}

            # Compress raw BLAST hits (kept for inspection)
            gzip -f {params.outdir}/silva_blast_hits.tsv
            """

rule assign_asv_taxonomy:
    """DADA2 RDP taxonomy assignment for ASVs, and the final ASV tables"""
    input:
        seqtab=get_asv_taxonomy_input,
    output:
        taxa=OUTPUT_DIR + "/03_dada2/taxID.rds",
        taxa_boot=OUTPUT_DIR + "/03_dada2/taxID_bootstrap.rds",
        combined_taxa=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        combined_taxa_boot=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_bootstrap_{DB_NAME}.txt",
    params:
        outdir=OUTPUT_DIR + "/03_dada2",
        amplicon=AMPLICON,
        taxonomy_db=EFFECTIVE_TAXONOMY_DB,
        project_name=PROJECT_NAME,
        db_name=DB_NAME,
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/assign_asv_taxonomy.log",
    threads: workflow.cores
    shell:
        """
        Rscript {params.pipeline_dir}/workflow/scripts/assign_asv_taxonomy.R \
          {input.seqtab} \
          {params.outdir} \
          {params.amplicon} \
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
        project_name=PROJECT_NAME,
        db_name=DB_NAME,
    shell:
        """
        cat > {output} << 'EOF'
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
- `{params.project_name}__ASV_counts_unfiltered.txt`: Counts for all ASVs after chimera removal, no taxonomy (input to the OTU pipeline)

**AMF SILVA screen (18S-AMF only, `amf_filter/`):**
- ASVs are BLASTed against SILVA before taxonomy; only those matching AMF reach the combined results and taxonomy files above
- `seqtab_nochim.rds`, `sequences.fasta` and the unfiltered counts keep all ASVs, and OTUs are clustered from all ASVs
- `amf_filter_asv_decisions.tsv`: Every ASV with its decision, reason, reads and best SILVA hit
- `amf_filter_discarded.tsv`: Discarded ASVs only
- `amf_filter_per_sample.tsv`: Reads kept / total per sample
- `amf_filter_counts.txt`: Summary counts (also reported in qc_summary.txt)
- `ASVs_AMF_filtered.fasta`: Kept ASV sequences
- `seqtab_nochim_AMF_filtered.rds`: Sequence table of kept ASVs
- `silva_blast_hits.tsv.gz`: Raw BLAST hits

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
        combined_taxa = OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
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
        pipeline_dir = workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/qc.log",
    shell:
        """
        mkdir -p {params.outdir}/figures
        Rscript {params.pipeline_dir}/workflow/scripts/qc_figures.R \
          {params.outdir} {params.amplicon} {params.platform} \
          {input.dada2_summary} {params.dada2_dir} \
          {params.primer_trimmed_dir} \
          2>&1 | tee {log}
        """
