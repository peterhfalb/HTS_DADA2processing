"""
Cutadapt rules: adapter and primer trimming
"""

# ============================================================================
# Step 1a: Adapter Trimming
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
    shell:
        """
        mkdir -p $(dirname {output.log})
        cutadapt --cores 8 --pair-filter=any --minimum-length 150 \
          -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA \
          -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT \
          -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT \
          --json={output.json} \
          -o {output.r1} \
          -p {output.r2} \
          {input.r1} {input.r2} \
          > {output.log} 2>&1
        """

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
    log:
        OUTPUT_DIR + "/.logs/primer_summary.log",
    shell:
        """
        grep "passing" {input.logs} > {output.summary} 2>&1 || true
        cat > {output.readme} << 'EOF'
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
- Amplicon type: {amplicon}

### Next step:
DADA2 denoising and taxonomy assignment in 03_dada2/
EOF
        """
