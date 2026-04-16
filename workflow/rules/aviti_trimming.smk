"""
Aviti trimming rules: 3-step Trimmomatic-cutadapt-Trimmomatic workflow
"""

# ============================================================================
# Step 1: Trimmomatic - First Nextera adapter pass
# ============================================================================

rule aviti_trim_adapters_pass1:
    """Remove Aviti/Nextera adapters from raw paired-end reads (first pass)"""
    input:
        r1=get_input_r1,
        r2=get_input_r2,
    output:
        r1=OUTPUT_DIR + "/01_adapter/{sample}_R1_001.paired.fastq.gz",
        r2=OUTPUT_DIR + "/01_adapter/{sample}_R2_001.paired.fastq.gz",
        r1_unpaired=temp(OUTPUT_DIR + "/01_adapter/{sample}_R1_001.unpaired.fastq.gz"),
        r2_unpaired=temp(OUTPUT_DIR + "/01_adapter/{sample}_R2_001.unpaired.fastq.gz"),
        log=OUTPUT_DIR + "/01_adapter/01_logs/{sample}_trimmomatic_pass1.log.txt",
    params:
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/aviti_trimmomatic_pass1_{sample}.log",
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.log})
        ADAPTER_FILE="{params.pipeline_dir}/workflow/adapters/aviti_adapters.fa"

        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          ILLUMINACLIP:$ADAPTER_FILE:2:30:10:2:True \
          LEADING:3 TRAILING:3 MINLEN:100 \
          -trimlog {output.log} \
          2>&1 | tee {log}
        """

# ============================================================================
# Step 2: Cutadapt - Primer Trimming
# ============================================================================

rule aviti_trim_primers:
    """Remove amplicon-specific primers from Trimmomatic-trimmed reads"""
    input:
        r1=OUTPUT_DIR + "/01_adapter/{sample}_R1_001.paired.fastq.gz",
        r2=OUTPUT_DIR + "/01_adapter/{sample}_R2_001.paired.fastq.gz",
    output:
        r1=OUTPUT_DIR + "/01b_primer_trimmed/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/01b_primer_trimmed/{sample}_R2_001.fastq.gz",
        log=OUTPUT_DIR + "/01b_primer_trimmed/02_logs/cutadapt.{sample}.log.txt",
        json=OUTPUT_DIR + "/01b_primer_trimmed/02_logs/cutadapt.{sample}.json",
    params:
        fwd_primer=EFFECTIVE_FWD,
        rev_primer=EFFECTIVE_REV,
    log:
        OUTPUT_DIR + "/.logs/aviti_cutadapt_primers_{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.log})
        cutadapt --cores 4 --pair-filter=any --minimum-length 100 --discard-untrimmed \
          -g {params.fwd_primer} \
          -G {params.rev_primer} \
          --json={output.json} \
          -o {output.r1} \
          -p {output.r2} \
          {input.r1} {input.r2} \
          > {output.log} 2>&1
        """

# ============================================================================
# Step 3: Trimmomatic - Second Nextera adapter pass
# ============================================================================

rule aviti_trim_adapters_pass2:
    """Remove remaining Nextera adapters (second pass)"""
    input:
        r1=OUTPUT_DIR + "/01b_primer_trimmed/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/01b_primer_trimmed/{sample}_R2_001.fastq.gz",
    output:
        r1=OUTPUT_DIR + "/02_primer_trimmed/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/02_primer_trimmed/{sample}_R2_001.fastq.gz",
        r1_unpaired=temp(OUTPUT_DIR + "/02_primer_trimmed/{sample}_R1_001.unpaired.fastq.gz"),
        r2_unpaired=temp(OUTPUT_DIR + "/02_primer_trimmed/{sample}_R2_001.unpaired.fastq.gz"),
        log=OUTPUT_DIR + "/02_primer_trimmed/02_logs/{sample}_trimmomatic_pass2.log.txt",
    params:
        pipeline_dir=workflow.basedir,
    log:
        OUTPUT_DIR + "/.logs/aviti_trimmomatic_pass2_{sample}.log",
    threads: 8
    shell:
        """
        mkdir -p $(dirname {output.log})
        ADAPTER_FILE="{params.pipeline_dir}/workflow/adapters/aviti_adapters.fa"

        trimmomatic PE -threads {threads} \
          {input.r1} {input.r2} \
          {output.r1} {output.r1_unpaired} \
          {output.r2} {output.r2_unpaired} \
          ILLUMINACLIP:$ADAPTER_FILE:2:30:10:2:True \
          LEADING:3 TRAILING:3 MINLEN:100 \
          -trimlog {output.log} \
          2>&1 | tee {log}
        """

# ============================================================================
# Step 4: Summarize Trimmomatic Results
# ============================================================================

rule aviti_summarize_trimming:
    """Aggregate Trimmomatic and cutadapt trimming statistics"""
    input:
        trimmo1_stderr=expand(
            OUTPUT_DIR + "/.logs/aviti_trimmomatic_pass1_{sample}.log",
            sample=SAMPLES,
        ),
        cutadapt_logs=expand(
            OUTPUT_DIR + "/01b_primer_trimmed/02_logs/cutadapt.{sample}.log.txt",
            sample=SAMPLES,
        ),
        trimmo2_stderr=expand(
            OUTPUT_DIR + "/.logs/aviti_trimmomatic_pass2_{sample}.log",
            sample=SAMPLES,
        ),
    output:
        summary_pass1=OUTPUT_DIR + "/01_adapter/summary_trimmomatic_pass1.txt",
        summary_pass2=OUTPUT_DIR + "/02_primer_trimmed/summary_trimmomatic_pass2.txt",
        readme_adapter=OUTPUT_DIR + "/01_adapter/README.txt",
        readme_primer=OUTPUT_DIR + "/01b_primer_trimmed/README.txt",
        readme_final=OUTPUT_DIR + "/02_primer_trimmed/README.txt",
    shell:
        """
        # Summarize Trimmomatic pass 1 (grep stderr logs for summary stats)
        echo "# Trimmomatic Pass 1 (Adapter Removal)" > {output.summary_pass1}
        echo "Sample\tInput Pairs\tBoth Surviving\tForward Only\tReverse Only\tDropped" >> {output.summary_pass1}
        for log in {input.trimmo1_stderr}; do
            sample=$(basename "$log" .log | sed 's/aviti_trimmomatic_pass1_//')
            grep "Input Read Pairs:" "$log" | head -1 | sed "s/^/$sample\t/" >> {output.summary_pass1}
        done 2>/dev/null || true

        # Summarize Trimmomatic pass 2 (grep stderr logs for summary stats)
        echo "# Trimmomatic Pass 2 (Final Adapter Removal)" > {output.summary_pass2}
        echo "Sample\tInput Pairs\tBoth Surviving\tForward Only\tReverse Only\tDropped" >> {output.summary_pass2}
        for log in {input.trimmo2_stderr}; do
            sample=$(basename "$log" .log | sed 's/aviti_trimmomatic_pass2_//')
            grep "Input Read Pairs:" "$log" | head -1 | sed "s/^/$sample\t/" >> {output.summary_pass2}
        done 2>/dev/null || true

        # Write README for 01_adapter
        cat > {output.readme_adapter} << 'EOF'
# Adapter Trimming Results (Trimmomatic Pass 1)
## Directory: 01_adapter/

This directory contains paired-end reads after removal of Aviti Nextera adapters (first pass).

### Contents:
- `<sample>_R1_001.paired.fastq.gz`: Forward reads with adapters removed
- `<sample>_R2_001.paired.fastq.gz`: Reverse reads with adapters removed
- `01_logs/`: Trimmomatic trimming logs
- `summary_trimmomatic_pass1.txt`: Aggregate of read counts per sample

### Parameters:
- ILLUMINACLIP: Aviti Nextera adapters (seed=2, palindrome=30, simple=10)
- keepBothReads: True (preserve both reads during palindrome detection)
- LEADING:3, TRAILING:3: Quality trimming from ends
- MINLEN:100: Minimum read length

### Next step:
Primer trimming in 01b_primer_trimmed/
EOF

        # Write README for 01b_primer_trimmed
        cat > {output.readme_primer} << 'EOF'
# Primer Trimming Results (Cutadapt)
## Directory: 01b_primer_trimmed/

This directory contains reads after removal of amplicon-specific primers.

### Contents:
- `<sample>_R1_001.fastq.gz`: Forward reads with primers removed
- `<sample>_R2_001.fastq.gz`: Reverse reads with primers removed
- `02_logs/`: Cutadapt trimming logs (json and text)

### Parameters:
- Minimum length: 100 bp
- Discard untrimmed: Yes (require both primers to be found)

### Next step:
Final adapter removal in 02_primer_trimmed/
EOF

        # Write README for 02_primer_trimmed
        cat > {output.readme_final} << 'EOF'
# Final Trimming Results (Trimmomatic Pass 2)
## Directory: 02_primer_trimmed/

This directory contains fully trimmed reads ready for DADA2 denoising.

### Contents:
- `<sample>_R1_001.fastq.gz`: Forward reads (final quality trimmed)
- `<sample>_R2_001.fastq.gz`: Reverse reads (final quality trimmed)
- `02_logs/`: Trimmomatic pass 2 logs
- `summary_trimmomatic_pass2.txt`: Final read counts per sample

### Parameters:
- ILLUMINACLIP: Aviti Nextera adapters (final pass)
- LEADING:3, TRAILING:3: Final quality trimming
- MINLEN:100: Minimum read length

### Next step:
DADA2 denoising pipeline in 03_dada2/
EOF
        """
