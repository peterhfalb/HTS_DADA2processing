"""
Cutadapt rules: adapter and primer trimming (Illumina platform)
"""

# ============================================================================
# Step 1a: Adapter Trimming
# ============================================================================

rule trim_adapters:
    """Remove Illumina adapters from paired-end reads"""
    input:
        r1=get_input_r1,
        r2=get_input_r2,
    output:
        r1=OUTPUT_DIR + "/01_adapter/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/01_adapter/{sample}_R2_001.fastq.gz",
        log=OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.log.txt",
        json=OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.json",
    log:
        OUTPUT_DIR + "/.logs/cutadapt_adapter_{sample}.log",
    threads: 8
    run:
        # Build adapter flags from detected adapters
        adapter_flags = ""
        adapters = DETECTED_ADAPTERS["adapters"]
        # Apply each detected adapter to both R1 (-a/-g) and R2 (-A/-G)
        # Assuming 3' end adapters, add to both -a and -A flags
        for adapter in adapters:
            adapter_flags += f" -a {adapter} -A {adapter}"

        shell(f"""
        mkdir -p $(dirname {output.log})
        cutadapt --cores {threads} --pair-filter=any --minimum-length 150 {adapter_flags} \
          --json={output.json} \
          -o {output.r1} \
          -p {output.r2} \
          {input.r1} {input.r2} \
          > {output.log} 2>&1
        """)

rule summarize_adapter_logs:
    """Aggregate adapter trimming log statistics"""
    input:
        logs=expand(
            OUTPUT_DIR + "/01_adapter/01_logs/cutadapt.{sample}.log.txt",
            sample=SAMPLES,
        ),
    output:
        summary=OUTPUT_DIR + "/01_adapter/summary_adapter_trimming.txt",
        readme=OUTPUT_DIR + "/01_adapter/README.txt",
    log:
        OUTPUT_DIR + "/.logs/adapter_summary.log",
    shell:
        """
        grep "passing" {input.logs} > {output.summary} 2>&1 || true
        cat > {output.readme} << 'EOF'
# Adapter Trimming Results
## Directory: 01_adapter/

This directory contains paired-end reads after removal of Illumina adapters.

### Contents:
- `<sample>_R1_001.fastq.gz`: Forward reads with adapters removed
- `<sample>_R2_001.fastq.gz`: Reverse reads with adapters removed
- `01_logs/`: Individual cutadapt logs
- `summary_adapter_trimming.txt`: Aggregate of "passing" reads per sample

### Parameters:
- Min length: 150 bp
- Adapters removed from 3' ends and reverse complements
- Discarded if either read pair falls below min length

### Next step:
Primer trimming in 02_primer_trimmed/
EOF
        """

# ============================================================================
# Step 1b: Primer Trimming
# ============================================================================

rule trim_primers:
    """Remove amplicon-specific primers"""
    input:
        r1=OUTPUT_DIR + "/01_adapter/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/01_adapter/{sample}_R2_001.fastq.gz",
    output:
        r1=OUTPUT_DIR + "/02_primer_trimmed/{sample}_R1_001.fastq.gz",
        r2=OUTPUT_DIR + "/02_primer_trimmed/{sample}_R2_001.fastq.gz",
        log=OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.log.txt",
        json=OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.json",
    params:
        fwd_primer=EFFECTIVE_FWD,
        rev_primer=EFFECTIVE_REV,
    log:
        OUTPUT_DIR + "/.logs/cutadapt_primers_{sample}.log",
    threads: 4
    shell:
        """
        mkdir -p $(dirname {output.log})
        cutadapt --cores {threads} --pair-filter=any --minimum-length 100 --discard-untrimmed \
          -g {params.fwd_primer} \
          -G {params.rev_primer} \
          --json={output.json} \
          -o {output.r1} \
          -p {output.r2} \
          {input.r1} {input.r2} \
          > {output.log} 2>&1
        """

rule summarize_primer_logs:
    """Aggregate primer trimming log statistics"""
    input:
        logs=expand(
            OUTPUT_DIR + "/02_primer_trimmed/02_logs/cutadapt.{sample}.log.txt",
            sample=SAMPLES,
        ),
    output:
        summary=OUTPUT_DIR + "/02_primer_trimmed/summary_primer_trimming.txt",
        readme=OUTPUT_DIR + "/02_primer_trimmed/README.txt",
    params:
        amplicon=AMPLICON,
    log:
        OUTPUT_DIR + "/.logs/primer_summary.log",
    shell:
        """
        grep "passing" {input.logs} > {output.summary} 2>&1 || true
        cat > {output.readme} << EOF
# Primer Trimming Results
## Directory: 02_primer_trimmed/

This directory contains paired-end reads after removal of amplicon-specific primers.
Reads without both primers have been discarded.

### Contents:
- `<sample>_R1_001.fastq.gz`: Forward reads with primers removed
- `<sample>_R2_001.fastq.gz`: Reverse reads with primers removed
- `02_logs/`: Individual cutadapt logs
- `summary_primer_trimming.txt`: Aggregate of "passing" reads per sample

### Parameters:
- Min length: 100 bp
- Discard untrimmed: Yes (requires both primers present)
- Amplicon type: {params.amplicon}

### Next step:
DADA2 denoising and taxonomy assignment in 03_dada2/
EOF
        """
