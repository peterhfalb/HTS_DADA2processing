"""
ASV to OTU Pipeline: VSEARCH clustering + mumu curation + taxonomy assignment
"""

# Helper function to choose input for VSEARCH (ITSx conditional)
def get_vsearch_input(wildcards):
    if RUN_ITSX and AMPLICON in ("ITS1", "ITS2"):
        return OUTPUT_DIR + f"/05_asv2otu/02_itsx/Centroid.ITSx.{AMPLICON}.filtered.fasta"
    return OUTPUT_DIR + "/05_asv2otu/01_input/Centroid.fasta"

# ============================================================================
# Step 1: Prepare OTU input from ASV table
# ============================================================================

rule prepare_otu_input:
    """Convert ASV table to FASTA and abundance table for OTU clustering"""
    input:
        asv_table=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
    output:
        fasta=OUTPUT_DIR + "/05_asv2otu/01_input/Centroid.fasta",
        abund=OUTPUT_DIR + "/05_asv2otu/01_input/asv_abundance_table.txt",
        map=OUTPUT_DIR + "/05_asv2otu/01_input/Map_file.csv",
    params:
        outdir=OUTPUT_DIR + "/05_asv2otu/01_input",
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/prepare_otu_input.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.pipeline_dir}/workflow/scripts/prepare_otu_input.R {input.asv_table} {params.outdir} 2>&1 | tee {log}
        """

# ============================================================================
# Step 2: ITSx extraction (ITS1/ITS2 only, optional)
# ============================================================================

# Only define itsx_extract rule if RUN_ITSX is True and we have an ITS amplicon
if RUN_ITSX and AMPLICON in ("ITS1", "ITS2"):
    rule itsx_extract:
        input:
            fasta=OUTPUT_DIR + "/05_asv2otu/01_input/Centroid.fasta",
        output:
            fasta=OUTPUT_DIR + f"/05_asv2otu/02_itsx/Centroid.ITSx.{AMPLICON}.filtered.fasta",
        params:
            outdir=OUTPUT_DIR + "/05_asv2otu/02_itsx",
            region=AMPLICON,
        log:
            OUTPUT_DIR + "/.logs/itsx.log",
        threads: workflow.cores
        shell:
            """
            mkdir -p {params.outdir}

            # Extract ITS region
            ITSx -i {input.fasta} \
                --complement F \
                -t F \
                --preserve T \
                --cpu {threads} \
                -o {params.outdir}/Centroid.ITSx 2>&1 | tee {log}

            # Remove zero-abundance sequences (filter out sequences with size=0 annotations)
            awk '/^>/ {{
                if (match($0, /;size=([0-9]+)/, a)) {{
                    size = a[1]
                }} else {{
                    size = 1  # Default to keep if no size annotation
                }}
                print_seq = (size > 0)
            }} !/^>/ {{
                if (print_seq) print
            }}' {params.outdir}/Centroid.ITSx.{params.region}.fasta > {params.outdir}/Centroid.ITSx.{params.region}.filtered.fasta

            # Clean up intermediate files (keep only the final filtered output)
            rm -f {params.outdir}/Centroid.ITSx.{params.region}.fasta
            """

# ============================================================================
# Step 3: VSEARCH clustering and chimera removal
# ============================================================================

rule vsearch_cluster:
    """VSEARCH: sort, cluster, chimera check, and build OTU table"""
    input:
        fasta=get_vsearch_input,
    output:
        centroids=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/{PROJECT_NAME}.centroids",
        otutable=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/{PROJECT_NAME}.otutable",
        stats=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/otu_processing_summary.txt",
    params:
        outdir=OUTPUT_DIR + "/05_asv2otu/03_vsearch",
        cluster_id=OTU_CLUSTER_ID,
        project=PROJECT_NAME,
        abund_table=OUTPUT_DIR + "/05_asv2otu/01_input/asv_abundance_table.txt",
        pipeline_dir=workflow.basedir,
    threads: workflow.cores
    log:
        OUTPUT_DIR + "/.logs/vsearch_cluster.log",
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        PIPELINE_DIR={params.pipeline_dir}

        # Count input sequences (already ITSx-filtered if RUN_ITSX=true)
        INPUT_ASVS=$(grep -c "^>" {input.fasta})

        # Step 1: Sort by abundance (size-based)
        vsearch --fasta_width 0 \
            --sortbysize {input.fasta} \
            --minseqlength 10 \
            --sizein --sizeout \
            --threads {threads} \
            --output sorted.fasta 2>&1 | tee cluster.log

        # Step 2: Cluster at specified identity
        vsearch --cluster_size sorted.fasta \
            --id {params.cluster_id} \
            --sizein --sizeout \
            --minseqlength 10 \
            --qmask none \
            --threads {threads} \
            --centroids clustered.fasta 2>&1 | tee -a cluster.log

        # Count post-cluster
        POST_CLUSTER=$(grep -c "^>" clustered.fasta)

        # Step 3: De novo chimera checking
        vsearch --uchime_denovo clustered.fasta \
            --sizein --sizeout \
            --qmask none \
            --threads {threads} \
            --nonchimeras nochim.fasta 2>&1 | tee -a cluster.log

        # Count post-chimera
        POST_CHIMERA=$(grep -c "^>" nochim.fasta)

        # Step 4: Final sort and save centroids
        vsearch --fasta_width 0 \
            --sortbysize nochim.fasta \
            --sizein --sizeout \
            --minseqlength 10 \
            --threads {threads} \
            --output {params.project}.centroids 2>&1 | tee -a cluster.log

        # Step 5: Map original sequences back to centroids
        vsearch --usearch_global {input.fasta} \
            --db {params.project}.centroids \
            --minseqlength 10 \
            --strand plus \
            --id {params.cluster_id} \
            --maxaccepts 0 \
            --qmask none \
            --threads {threads} \
            --uc {params.project}.uc 2>&1 | tee -a cluster.log

        # Step 6: Build OTU table using R script
        Rscript $PIPELINE_DIR/workflow/scripts/build_otu_table.R \
            {params.project}.uc \
            {params.abund_table} \
            {params.project}.otutable

        # Calculate and save OTU clustering statistics
        CHIMERAS_REMOVED=$((POST_CLUSTER - POST_CHIMERA))
        PCT_CHIMERIC=$(awk -v cr=$CHIMERAS_REMOVED -v pc=$POST_CLUSTER 'BEGIN {{printf "%.1f", (cr / pc) * 100}}')
        PCT_ASVS_REMOVED=$(awk -v ia=$INPUT_ASVS -v pc=$POST_CLUSTER 'BEGIN {{printf "%.1f", ((ia - pc) / ia) * 100}}')

        # Save all counts to file before cleanup (needed by OTU QC summary script)
        echo "INPUT_ASVS=$INPUT_ASVS" > {params.project}_vsearch_counts.txt
        echo "POST_CLUSTER=$POST_CLUSTER" >> {params.project}_vsearch_counts.txt
        echo "CHIMERAS_REMOVED=$CHIMERAS_REMOVED" >> {params.project}_vsearch_counts.txt
        echo "POST_CHIMERA=$POST_CHIMERA" >> {params.project}_vsearch_counts.txt
        echo "PCT_CHIMERIC=$PCT_CHIMERIC" >> {params.project}_vsearch_counts.txt
        echo "PCT_ASVS_REMOVED=$PCT_ASVS_REMOVED" >> {params.project}_vsearch_counts.txt

        # Clean up intermediate files
        rm -f sorted.fasta clustered.fasta nochim.fasta {params.project}.uc

        # Copy counts file to output
        cp {params.project}_vsearch_counts.txt {output.stats}

        echo "VSEARCH clustering complete" >> cluster.log
        """

# ============================================================================
# Step 4: mumu curation
# ============================================================================

rule mumu_curation:
    """mumu: OTU curation based on read abundance and BLAST matches"""
    input:
        centroids=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/{PROJECT_NAME}.centroids",
        otutable=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/{PROJECT_NAME}.otutable",
    output:
        curated_table=OUTPUT_DIR + f"/05_asv2otu/04_mumu/{PROJECT_NAME}_mumu_curated.txt",
        curated_fasta=OUTPUT_DIR + f"/05_asv2otu/04_mumu/Centroid_mumu_curated.fas",
    params:
        outdir=OUTPUT_DIR + "/05_asv2otu/04_mumu",
        project=PROJECT_NAME,
        blast_id=MUMU_BLAST_ID,
        mumu_ratio=MUMU_RATIO,
    threads: workflow.cores
    log:
        OUTPUT_DIR + "/.logs/mumu_curation.log",
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}

        # Verify mumu binary exists
        if [[ ! -f "$HOME/packages/mumu/mumu" ]]; then
            echo "ERROR: mumu binary not found at $HOME/packages/mumu/mumu"
            echo "Please run install.sh to compile mumu"
            exit 1
        fi

        # Strip size= annotations from centroids for BLAST
        sed 's/;size=[0-9]*//g' {input.centroids} > OTU_centroids

        # Create BLAST database
        makeblastdb -in OTU_centroids -parse_seqids -dbtype nucl 2>&1 | tee {log}

        # Run BLAST all-vs-all with mumu parameters
        blastn -db OTU_centroids \
            -outfmt '6 qseqid sseqid pident' \
            -out match_list.txt \
            -qcov_hsp_perc 80 \
            -perc_identity {params.blast_id} \
            -num_threads {threads} \
            -query OTU_centroids 2>&1 | tee -a {log}

        # Prepare OTU table for mumu (strip ;size= from OTU_ID, preserve header)
        awk 'BEGIN{{FS=OFS="\\t"}} {{sub(/;.*/, "", $1); print}}' {input.otutable} > mumu_table.txt

        # Run mumu curation
        $HOME/packages/mumu/mumu \
            --otu_table mumu_table.txt \
            --match_list match_list.txt \
            --log mumu.log \
            --new_otu_table {params.project}_mumu_curated.txt \
            --minimum_ratio {params.mumu_ratio} 2>&1 | tee -a {log}

        # Extract curated centroid sequences from original centroids
        # Get list of OTU IDs that passed mumu curation
        awk 'NR > 1 {{print $1}}' {params.project}_mumu_curated.txt > otu_ids.txt

        # For each OTU, extract its sequence and calculate total abundance for size annotation
        > Centroid_mumu_curated.fas
        while IFS= read -r otu_id; do
            # Get sequence header and sequence from original centroids
            grep "^>$otu_id;size" {input.centroids} -A1 | head -2 >> Centroid_mumu_curated.fas
        done < otu_ids.txt

        # Clean up intermediate files (output files are already in the output directory)
        rm -f OTU_centroids* match_list.txt otu_ids.txt mumu_table.txt

        echo "mumu curation complete" >> {log}
        """

# ============================================================================
# Step 5: Assign OTU taxonomy
# ============================================================================

rule assign_otu_taxonomy:
    """DADA2 RDP taxonomy assignment for OTU centroids"""
    input:
        fasta=OUTPUT_DIR + f"/05_asv2otu/04_mumu/Centroid_mumu_curated.fas",
    output:
        taxonomy=OUTPUT_DIR + f"/05_asv2otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}_combined.txt",
    params:
        outdir=OUTPUT_DIR + "/05_asv2otu/05_taxonomy",
        taxonomy_db=EFFECTIVE_TAXONOMY_DB,
        amplicon=AMPLICON,
        db_name=DB_NAME,
        pipeline_dir=workflow.basedir,
        output_base=OUTPUT_DIR + f"/05_asv2otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}",
    threads: workflow.cores
    log:
        OUTPUT_DIR + "/.logs/assign_otu_taxonomy.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.pipeline_dir}/workflow/scripts/assign_otu_taxonomy.R \
            {input.fasta} \
            {params.taxonomy_db} \
            {params.output_base}.txt \
            {params.amplicon} \
            {params.db_name} \
            2>&1 | tee {log}
        """

# ============================================================================
# Step 6: Combine OTU table with taxonomy
# ============================================================================

rule combine_otu_taxonomy:
    """Merge OTU abundances with taxonomy assignments"""
    input:
        otu_table=OUTPUT_DIR + f"/05_asv2otu/04_mumu/{PROJECT_NAME}_mumu_curated.txt",
        taxonomy=OUTPUT_DIR + f"/05_asv2otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}_combined.txt",
        map_file=OUTPUT_DIR + "/05_asv2otu/01_input/Map_file.csv",
    output:
        combined=OUTPUT_DIR + f"/05_asv2otu/{PROJECT_NAME}__OTUs_with_taxonomy_{DB_NAME}.txt",
    params:
        outdir=OUTPUT_DIR + "/05_asv2otu",
        project=PROJECT_NAME,
        amplicon=AMPLICON,
        db_name=DB_NAME,
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/combine_otu_taxonomy.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.pipeline_dir}/workflow/scripts/combine_otu_taxonomy.R \
            {params.project} \
            {params.amplicon} \
            {input.otu_table} \
            {input.taxonomy} \
            {input.map_file} \
            {output.combined} \
            2>&1 | tee {log}
        """

# ============================================================================
# Step 7: Generate OTU QC summary
# ============================================================================

rule generate_otu_qc_summary:
    """Generate comprehensive OTU QC summary from all pipeline steps"""
    input:
        asv_table=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        centroids=OUTPUT_DIR + f"/05_asv2otu/03_vsearch/{PROJECT_NAME}.centroids",
        mumu_table=OUTPUT_DIR + f"/05_asv2otu/04_mumu/{PROJECT_NAME}_mumu_curated.txt",
        taxonomy=OUTPUT_DIR + f"/05_asv2otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}_combined.txt",
    output:
        qc_summary=OUTPUT_DIR + "/05_asv2otu/otu_qc_summary.txt",
    params:
        output_dir=OUTPUT_DIR,
        project=PROJECT_NAME,
        run_itsx="1" if (RUN_ITSX and AMPLICON in ("ITS1", "ITS2")) else "0",
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/generate_otu_qc_summary.log",
    shell:
        """
        Rscript {params.pipeline_dir}/workflow/scripts/generate_otu_qc_summary.R \
            {params.output_dir} \
            {params.project} \
            {params.run_itsx} \
            2>&1 | tee {log}
        """

# ============================================================================
# Step 7b: Generate READMEs for each asv2otu subfolder
# ============================================================================

rule asv2otu_input_readme:
    """Create README for OTU input preparation step"""
    output:
        OUTPUT_DIR + "/05_asv2otu/01_input/README.txt",
    shell:
        """
        cat > {output} << 'EOF'
# ASV to OTU Pipeline - Step 1: Input Preparation
## Directory: 05_asv2otu/01_input/

This directory contains the prepared input files for OTU clustering, derived from the DADA2 ASV table.

### Contents:

**Input Files:**
- `Centroid.fasta`: FASTA sequences of all ASVs with size annotations (;size=N)
- `asv_abundance_table.txt`: Tab-separated abundance table with OTU IDs and sample counts
- `Map_file.csv`: Sample metadata mapping file

### Processing Notes:
- All sequences retain abundance annotations for VSEARCH clustering
- OTU IDs are derived from ASV sequences (Seq001, Seq002, etc.)
- Sequences are sorted by decreasing abundance for optimal clustering

### Next Step:
Sequences will be sorted by abundance, clustered at 97% identity, and screened for chimeras using VSEARCH.
EOF
        """

rule asv2otu_itsx_readme:
    """Create README for ITSx extraction step (ITS amplicons only)"""
    output:
        OUTPUT_DIR + "/05_asv2otu/02_itsx/README.txt",
    shell:
        """
        cat > {output} << 'EOF'
# ASV to OTU Pipeline - Step 2: ITSx Extraction
## Directory: 05_asv2otu/02_itsx/

This directory contains ITSx-extracted sequences for ITS amplicons (ITS1 or ITS2).

### Contents:

**Filtered Sequences:**
- `Centroid.ITSx.[REGION].filtered.fasta`: ITS region-specific sequences with size annotations

### Processing Notes:
- ITSx extracts conserved ITS regions, removing 18S and 28S rRNA regions
- Only sequences with identified ITS region are retained
- Size annotations are preserved for downstream clustering
- This step applies only to ITS1 and ITS2 amplicons; other amplicon types skip this step

### Next Step:
Filtered sequences proceed to VSEARCH clustering at 97% identity.
EOF
        """

rule asv2otu_vsearch_readme:
    """Create README for VSEARCH clustering step"""
    output:
        OUTPUT_DIR + "/05_asv2otu/03_vsearch/README.txt",
    params:
        project_name=PROJECT_NAME,
    shell:
        """
        cat > {output} << EOF
# ASV to OTU Pipeline - Step 3: VSEARCH Clustering & Chimera Removal
## Directory: 05_asv2otu/03_vsearch/

This directory contains OTU sequences after clustering and chimera detection.

### Contents:

**Output Files:**
- `{params.project_name}.centroids`: Representative sequences for each OTU (after clustering, before chimera removal)
- `{params.project_name}.otutable`: OTU abundance table with sample counts (post-chimera removal)
- `otu_processing_summary.txt`: QC statistics for clustering and chimera removal steps

**QC Metrics (from otu_processing_summary.txt):**
- Input ASVs: Number of unique sequences entering clustering
- OTUs after clustering (97% identity): Sequences clustered at 97% identity
- Chimeras removed: Count of chimeric sequences detected and removed
- OTUs after chimera removal: Final count before mumu curation
- % chimeric: Percentage of clustered OTUs that were detected as chimeras

### Processing Details:
- **Abundance sorting:** Sequences sorted by decreasing abundance (VSEARCH --sortbysize)
- **Clustering:** Greedy clustering at 97% identity (VSEARCH --cluster_size)
- **Chimera detection:** De novo chimera detection with size-aware scoring
- **Mapping:** All original ASVs mapped back to OTU centroids for abundance accounting

### Next Step:
OTU centroids and abundance table proceed to mumu curation for redundancy reduction.
EOF
        """

rule asv2otu_mumu_readme:
    """Create README for mumu curation step"""
    output:
        OUTPUT_DIR + "/05_asv2otu/04_mumu/README.txt",
    params:
        project_name=PROJECT_NAME,
    shell:
        """
        cat > {output} << EOF
# ASV to OTU Pipeline - Step 4: mumu Curation
## Directory: 05_asv2otu/04_mumu/

This directory contains OTU sequences and abundance table after mumu-based curation and redundancy reduction.

### Contents:

**Output Files:**
- `{params.project_name}_mumu_curated.txt`: Curated OTU abundance table with mumu-selected OTUs only
- `Centroid_mumu_curated.fas`: Representative sequences for curated OTUs with abundance annotations

### Processing Details:
- **Redundancy reduction:** mumu identifies and merges OTUs that are very similar
  - BLAST all-vs-all comparison at specified identity (default 84% for non-16S-V4, 94% for 16S-V4)
  - Merges lower-abundance OTUs into higher-abundance ones within similarity threshold
  - Applies minimum abundance ratio threshold to prevent over-merging

- **Quality control:** Retains only OTUs passing curation filters
  - Removes singleton OTUs and low-abundance artifacts
  - Consolidates similar sequences to improve ecological interpretation

### Next Step:
Curated OTU sequences proceed to RDP taxonomy assignment.
EOF
        """

rule asv2otu_taxonomy_readme:
    """Create README for taxonomy assignment step"""
    output:
        OUTPUT_DIR + "/05_asv2otu/05_taxonomy/README.txt",
    shell:
        """
        cat > {output} << 'EOF'
# ASV to OTU Pipeline - Step 5: Taxonomy Assignment
## Directory: 05_asv2otu/05_taxonomy/

This directory contains taxonomy assignments for curated OTU centroids.

### Contents:

**Output Files:**
- `Taxonomy_rdp_[AMPLICON]_combined.txt`: Taxonomy assignments with bootstrap confidence values
  Columns: OTU_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species
           (plus corresponding bootstrap columns for confidence)

### Processing Details:
- **Database:** RDP (Ribosomal Database Project) or alternative taxonomy database
- **Method:** DADA2 assignTaxonomy with 100bp k-mer window
- **Bootstrap:** 100 iterations to estimate confidence for each taxonomy rank
- **Output:** Tab-separated format with full 7-rank taxonomy and bootstrap values

### Next Step:
Taxonomy assignments merge with curated OTU abundance table for final output.
EOF
        """

rule asv2otu_overview_readme:
    """Create overview README for entire ASV-to-OTU pipeline"""
    output:
        OUTPUT_DIR + "/05_asv2otu/README.txt",
    shell:
        """
        cat > {output} << 'EOF'
# ASV to OTU Pipeline Overview
## Directory: 05_asv2otu/

This directory contains the complete ASV-to-OTU clustering and curation workflow.

### Pipeline Steps:

**Step 1: Input Preparation** (01_input/)
- Converts DADA2 ASV table to FASTA and abundance formats
- Outputs: Centroid sequences, abundance table, sample map

**Step 2: ITSx Extraction** (02_itsx/)
- *ITS amplicons only:* Extracts conserved ITS1 or ITS2 regions
- Removes non-ITS sequences (18S, 28S rRNA regions)
- Outputs: Region-specific filtered sequences

**Step 3: VSEARCH Clustering & Chimera Removal** (03_vsearch/)
- Sorts sequences by abundance
- Clusters at 97% identity (operational taxonomic units)
- De novo chimera detection removes likely artifacts
- Outputs: Centroid sequences, OTU abundance table, QC statistics

**Step 4: mumu Curation** (04_mumu/)
- Reduces redundancy through abundance-weighted clustering
- BLAST-based matching identifies similar OTUs
- Merges low-abundance into high-abundance OTUs
- Outputs: Curated OTU abundance table and representative sequences

**Step 5: Taxonomy Assignment** (05_taxonomy/)
- RDP-based taxonomy assignment for OTU centroids
- Bootstrap confidence values for each taxonomic rank
- Outputs: OTU taxonomy table with confidence scores

### Final Outputs:

- **[PROJECT]__OTUs_with_taxonomy_[DB].txt**: Combined OTU abundance table with full taxonomy
- **otu_qc_summary.txt**: Quality control metrics from all clustering and curation steps

### Key Metrics Tracked:

- Input ASVs entering pipeline
- OTUs after 97% clustering
- Chimeras removed (count and percentage)
- OTUs removed by mumu curation
- OTUs with successful taxonomy assignment
- Processing percentages at each step

### Configuration Parameters:

- `--otu-cluster-id`: Clustering identity threshold (default 0.97 for 97%)
- `--mumu-blast-id`: BLAST identity for mumu matching (default varies by amplicon)
- `--mumu-ratio`: Minimum abundance ratio for mumu merging (default varies by amplicon)
- `--run-itsx`: Enable ITSx extraction for ITS amplicons (default: disabled)

### Interpretation:

- **OTU centroids:** Representative sequences at 97% identity
- **mumu curation:** Reduces false OTUs while preserving true diversity
- **Taxonomy:** Provides functional classification for downstream analysis
- **QC metrics:** Allow assessment of sequence quality and processing efficiency

### Next Analysis:

Use the combined OTU abundance table for:
- Alpha and beta diversity analysis
- Comparative abundance across samples
- Phylogenetic analysis and tree building
- Statistical testing and visualization
EOF
        """

# ============================================================================
# Step 8: Finalize outputs (copy to main dir, combine QC tables, write README)
# ============================================================================

rule finalize_outputs:
    """Copy summary files to main directory and generate final README"""
    input:
        asv_table=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_bootstrap_{DB_NAME}.txt",
        qc_per_sample=OUTPUT_DIR + "/04_dada2_QCsummary/qc_summary.txt",
        # ASV2OTU README files (conditional on SKIP_OTU)
        **{f"readme_input": OUTPUT_DIR + "/05_asv2otu/01_input/README.txt"} if not SKIP_OTU else {},
        **{f"readme_itsx": OUTPUT_DIR + "/05_asv2otu/02_itsx/README.txt"} if (not SKIP_OTU and RUN_ITSX and AMPLICON in ("ITS1", "ITS2")) else {},
        **{f"readme_vsearch": OUTPUT_DIR + "/05_asv2otu/03_vsearch/README.txt"} if not SKIP_OTU else {},
        **{f"readme_mumu": OUTPUT_DIR + "/05_asv2otu/04_mumu/README.txt"} if not SKIP_OTU else {},
        **{f"readme_taxonomy": OUTPUT_DIR + "/05_asv2otu/05_taxonomy/README.txt"} if not SKIP_OTU else {},
        **{f"readme_overview": OUTPUT_DIR + "/05_asv2otu/README.txt"} if not SKIP_OTU else {},
        **{f"otu_table": OUTPUT_DIR + f"/05_asv2otu/{PROJECT_NAME}__OTUs_with_taxonomy_{DB_NAME}.txt"} if not SKIP_OTU else {},
        **{f"otu_qc_summary": OUTPUT_DIR + f"/05_asv2otu/otu_qc_summary.txt"} if not SKIP_OTU else {},
    output:
        qc_main=OUTPUT_DIR + "/qc_summary.txt",
        readme=OUTPUT_DIR + "/README.txt",
    params:
        output_dir=OUTPUT_DIR,
        project=PROJECT_NAME,
        db_name=DB_NAME,
        skip_otu=SKIP_OTU,
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/finalize_outputs.log",
    shell:
        """
        Rscript {params.pipeline_dir}/workflow/scripts/finalize_outputs.R \
            {params.output_dir} \
            {params.project} \
            {params.db_name} \
            {params.skip_otu} \
            2>&1 | tee {log}
        """
