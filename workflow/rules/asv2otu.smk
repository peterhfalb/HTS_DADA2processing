"""
ASV to OTU Pipeline: VSEARCH clustering + mumu curation + taxonomy assignment
"""

# Helper function to choose input for VSEARCH (ITSx conditional)
def get_vsearch_input(wildcards):
    if RUN_ITSX and AMPLICON in ("ITS1", "ITS2"):
        return OUTPUT_DIR + f"/05_otu/02_itsx/Centroid.ITSx.{AMPLICON}.filtered.fasta"
    return OUTPUT_DIR + "/05_otu/01_input/Centroid.fasta"

# ============================================================================
# Step 1: Prepare OTU input from ASV table
# ============================================================================

rule prepare_otu_input:
    """Convert ASV table to FASTA and abundance table for OTU clustering"""
    input:
        asv_table=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
    output:
        fasta=OUTPUT_DIR + "/05_otu/01_input/Centroid.fasta",
        abund=OUTPUT_DIR + "/05_otu/01_input/asv_abundance_table.txt",
        map=OUTPUT_DIR + "/05_otu/01_input/Map_file.csv",
    params:
        outdir=OUTPUT_DIR + "/05_otu/01_input",
    log:
        OUTPUT_DIR + "/.logs/prepare_otu_input.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript workflow/scripts/prepare_otu_input.R {input.asv_table} {params.outdir} 2>&1 | tee {log}
        """

# ============================================================================
# Step 2: ITSx extraction (ITS1/ITS2 only, optional)
# ============================================================================

# Only define itsx_extract rule if RUN_ITSX is True and we have an ITS amplicon
if RUN_ITSX and AMPLICON in ("ITS1", "ITS2"):
    rule itsx_extract:
        input:
            fasta=OUTPUT_DIR + "/05_otu/01_input/Centroid.fasta",
        output:
            fasta=OUTPUT_DIR + f"/05_otu/02_itsx/Centroid.ITSx.{AMPLICON}.filtered.fasta",
        params:
            outdir=OUTPUT_DIR + "/05_otu/02_itsx",
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

            # Fix header format (ensure trailing semicolons)
            sed 's/;size=/;size=/g' {params.outdir}/Centroid.ITSx.{params.region}.fasta > {params.outdir}/Centroid.ITSx.{params.region}.fixed.fasta

            # Remove zero-abundance sequences
            awk '/^>/ {{
                if (match($0, /;size=([0-9]+)/, a)) {{
                    size = a[1]
                    if (size > 0) print
                    else getline
                }} else print
            }} !/^>/ {{
                if (size > 0) print
            }}' {params.outdir}/Centroid.ITSx.{params.region}.fixed.fasta > {params.outdir}/Centroid.ITSx.{params.region}.filtered.fasta

            rm -f {params.outdir}/Centroid.ITSx.* {params.outdir}/Centroid.ITSx.{params.region}.fixed.fasta
            """

# ============================================================================
# Step 3: VSEARCH clustering and chimera removal
# ============================================================================

rule vsearch_cluster:
    """VSEARCH: sort, cluster, chimera check, and build OTU table"""
    input:
        fasta=get_vsearch_input,
    output:
        centroids=OUTPUT_DIR + f"/05_otu/03_vsearch/{PROJECT_NAME}.centroids",
        otutable=OUTPUT_DIR + f"/05_otu/03_vsearch/{PROJECT_NAME}.otutable",
        stats=OUTPUT_DIR + f"/05_otu/03_vsearch/otu_processing_summary.txt",
    params:
        outdir=OUTPUT_DIR + "/05_otu/03_vsearch",
        cluster_id=OTU_CLUSTER_ID,
        project=PROJECT_NAME,
        abund_table=OUTPUT_DIR + "/05_otu/01_input/asv_abundance_table.txt",
    threads: workflow.cores
    log:
        OUTPUT_DIR + "/.logs/vsearch_cluster.log",
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}

        # Count input sequences
        INPUT_ASVS=$(grep -c "^>" {input.fasta})

        # Count post-ITSx (if ITSx was used)
        POST_ITSX="NA"
        if grep -q "ITSx" "{input.fasta}"; then
            POST_ITSX=$(grep -c "^>" {input.fasta})
        fi

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
        Rscript workflow/scripts/build_otu_table.R \
            {params.project}.uc \
            {params.abund_table} \
            {params.project}.otutable

        # Count classified OTUs (non-header lines in otutable)
        POST_MUMU=$(wc -l < {params.project}.otutable)
        POST_MUMU=$((POST_MUMU - 1))  # subtract header line

        # Write OTU processing summary
        echo -e "input_asvs\\tpost_itsx\\tpost_cluster\\tpost_chimera\\tpost_mumu\\tclassified_otus" > {params.project}_processing_summary.txt
        echo -e "$INPUT_ASVS\\t$POST_ITSX\\t$POST_CLUSTER\\t$POST_CHIMERA\\t$POST_MUMU\\t$POST_MUMU" >> {params.project}_processing_summary.txt

        # Clean up intermediate files
        rm -f sorted.fasta clustered.fasta nochim.fasta {params.project}.uc

        # Copy files to output
        cp {params.project}_processing_summary.txt {output.stats}

        echo "VSEARCH clustering complete" >> cluster.log
        """

# ============================================================================
# Step 4: mumu curation
# ============================================================================

rule mumu_curation:
    """mumu: OTU curation based on read abundance and BLAST matches"""
    input:
        centroids=OUTPUT_DIR + f"/05_otu/03_vsearch/{PROJECT_NAME}.centroids",
        otutable=OUTPUT_DIR + f"/05_otu/03_vsearch/{PROJECT_NAME}.otutable",
    output:
        curated_table=OUTPUT_DIR + f"/05_otu/04_mumu/{PROJECT_NAME}_mumu_curated.txt",
        curated_fasta=OUTPUT_DIR + f"/05_otu/04_mumu/Centroid_mumu_curated.fas",
    params:
        outdir=OUTPUT_DIR + "/05_otu/04_mumu",
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

        # Prepare OTU table for mumu (remove header, keep OTU_ID column)
        tail -n +2 {input.otutable} > mumu_table.txt

        # Run mumu curation
        $HOME/packages/mumu/mumu \
            --otu_table mumu_table.txt \
            --match_list match_list.txt \
            --log mumu.log \
            --new_otu_table {params.project}_mumu_curated.txt \
            --minimum_ratio {params.mumu_ratio} 2>&1 | tee -a {log}

        # Extract curated centroid sequences from original centroids
        # by grepping OTU IDs from curated table against original centroids
        awk 'NR > 1 {{print $1}}' {params.project}_mumu_curated.txt | \
            xargs -I {{}} grep "^>{{};size" {input.centroids} -A1 > Centroid_mumu_curated_temp.fas

        # Add size= annotations back from mumu curated table
        # This is complex, so we'll use a different approach:
        # Extract OTU IDs and their total abundances from curated table
        awk 'NR > 1 {{id=$1; abund=0; for(i=2; i<=NF; i++) abund+=$i; print id, abund}}' {params.project}_mumu_curated.txt > otu_abundances.txt

        # Build new FASTA with size annotations
        {{
            while IFS=' ' read otu_id abund; do
                grep "^>${{otu_id}};size" {input.centroids} | sed "s/;size=[0-9]*//;s/;size=/;size=${{abund}};/" | head -1
                grep "^>${{otu_id}};size" {input.centroids} -A1 | tail -1
            done < otu_abundances.txt
        }} > Centroid_mumu_curated.fas

        # Clean up
        rm -f OTU_centroids* match_list.txt otu_abundances.txt mumu_table.txt Centroid_mumu_curated_temp.fas

        # Copy output file
        cp {params.project}_mumu_curated.txt {output.curated_table}
        cp Centroid_mumu_curated.fas {output.curated_fasta}

        echo "mumu curation complete" >> {log}
        """

# ============================================================================
# Step 5: Assign OTU taxonomy
# ============================================================================

rule assign_otu_taxonomy:
    """DADA2 RDP taxonomy assignment for OTU centroids"""
    input:
        fasta=OUTPUT_DIR + f"/05_otu/04_mumu/Centroid_mumu_curated.fas",
    output:
        taxonomy=OUTPUT_DIR + f"/05_otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}_combined.txt",
    params:
        outdir=OUTPUT_DIR + "/05_otu/05_taxonomy",
        taxonomy_db=EFFECTIVE_TAXONOMY_DB,
        amplicon=AMPLICON,
        db_name=DB_NAME,
    threads: workflow.cores
    log:
        OUTPUT_DIR + "/.logs/assign_otu_taxonomy.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript workflow/scripts/assign_otu_taxonomy.R \
            {input.fasta} \
            {params.taxonomy_db} \
            {params.outdir} \
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
        otu_table=OUTPUT_DIR + f"/05_otu/04_mumu/{PROJECT_NAME}_mumu_curated.txt",
        taxonomy=OUTPUT_DIR + f"/05_otu/05_taxonomy/Taxonomy_rdp_{AMPLICON}_combined.txt",
        map_file=OUTPUT_DIR + "/05_otu/01_input/Map_file.csv",
    output:
        combined=OUTPUT_DIR + f"/05_otu/{PROJECT_NAME}__OTUs_with_taxonomy_{DB_NAME}.txt",
    params:
        outdir=OUTPUT_DIR + "/05_otu",
        project=PROJECT_NAME,
        amplicon=AMPLICON,
        db_name=DB_NAME,
    log:
        OUTPUT_DIR + "/.logs/combine_otu_taxonomy.log",
    shell:
        """
        mkdir -p {params.outdir}
        Rscript workflow/scripts/combine_otu_taxonomy.R \
            {params.project} \
            {params.amplicon} \
            {input.otu_table} \
            {input.taxonomy} \
            {input.map_file} \
            {output.combined} \
            2>&1 | tee {log}
        """

# ============================================================================
# Step 7: Finalize outputs (copy to main dir, combine QC tables, write README)
# ============================================================================

rule finalize_outputs:
    """Copy summary files to main directory and generate final README"""
    input:
        asv_table=OUTPUT_DIR + f"/03_dada2/{PROJECT_NAME}__combined_sequences_ASVtaxa_{DB_NAME}.txt",
        qc_per_sample=OUTPUT_DIR + "/04_QC/qc_summary.txt",
        **{f"otu_table": OUTPUT_DIR + f"/05_otu/{PROJECT_NAME}__OTUs_with_taxonomy_{DB_NAME}.txt"} if not SKIP_OTU else {},
        **{f"otu_stats": OUTPUT_DIR + f"/05_otu/03_vsearch/otu_processing_summary.txt"} if not SKIP_OTU else {},
    output:
        qc_main=OUTPUT_DIR + "/qc_summary.txt",
        readme=OUTPUT_DIR + "/README.txt",
    params:
        output_dir=OUTPUT_DIR,
        project=PROJECT_NAME,
        db_name=DB_NAME,
        skip_otu=SKIP_OTU,
    log:
        OUTPUT_DIR + "/.logs/finalize_outputs.log",
    shell:
        """
        Rscript workflow/scripts/finalize_outputs.R \
            {params.output_dir} \
            {params.project} \
            {params.db_name} \
            {params.skip_otu} \
            2>&1 | tee {log}
        """
