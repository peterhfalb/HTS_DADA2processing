# DADA2 Snakemake Pipeline

A flexible Snakemake-based pipeline for amplicon sequencing denoising and taxonomy assignment using DADA2. Designed for the Kennedy Lab at UMN, runs on the MSI Agate HPC cluster via SLURM.

## Supported Amplicon Types

- **16S-V4**: Bacterial and archaeal 16S rRNA (V4 region)
- **ITS1**: Fungal ITS1 region
- **ITS2**: Fungal ITS2 region
- **18S-AMF**: Arbuscular mycorrhizal fungi 18S
- **18S-V4**: Microbial eukaryotes 18S (V4 region)

## Quick Start

### 1. Installation

Add the pipeline to your PATH:
```bash
cd /path/to/DADA2_Snakemake_groundUp
bash install.sh
source ~/.bashrc
```

### 2. Run the Pipeline

```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir /path/to/raw/fastqs \
  --output-dir /path/to/output \
  --email your.email@umn.edu
```

Optional flags:
- `--quality good` or `--quality bad` (default: good)
- `--jobs N` for max parallel Snakemake jobs (default: 16)

### 3. Monitor the Job

The pipeline generates a SLURM submission script in your output directory. Submit with:
```bash
sbatch <output_dir>/submit_job.sh
```

Check logs:
```bash
tail -f <output_dir>/logs/snakemake_*.log
```

## Pipeline Overview

```
Input (raw fastqs)
        ↓
01_adapter         ← Illumina adapter trimming (8 cores)
        ↓
02_primer_trimmed  ← Amplicon-specific primer removal (4 cores)
        ↓
03_dada2           ← DADA2 denoising & taxonomy assignment
        ↓
Output (sequence tables, taxonomies, fasta)
```

### Step 1: Adapter Trimming
- **Tool**: cutadapt
- **Input**: Raw demultiplexed paired-end fastq files
- **Output**: `01_adapter/` with adapter-trimmed reads
- **Removes**: Universal Illumina adapters (3' ends)

### Step 2: Primer Trimming
- **Tool**: cutadapt
- **Input**: Adapter-trimmed reads
- **Output**: `02_primer_trimmed/` with primer-trimmed reads
- **Removes**: Amplicon-specific forward and reverse primers
- **Note**: Reads missing primers are discarded (`--discard-untrimmed`)

### Step 3: DADA2 Processing
- **Tool**: DADA2 (R package)
- **Steps**:
  1. Quality filtering (parameters depend on amplicon type and quality setting)
  2. Dereplication
  3. Error learning and DADA2 inference
  4. Paired-end merging
  5. Chimera removal (bimera, consensus method)
  6. Taxonomy assignment

**Output Directory: `03_dada2/`**
- `seqtab_nochim.rds`: Final sequence abundance table
- `sequences.fasta`: Unique sequences
- `<AMPLICON>_combined_sequences_taxa.txt`: Sequence counts with taxonomy
- `<AMPLICON>_combined_sequences_taxa_bootstrap.txt`: With bootstrap confidence
- `sequence_process_summary.txt`: Per-sample read counts at each step

## Input File Format

Input fastq files must follow the naming convention:
```
{ProjectName}_{SampleName}_L001_R1_001.fastq.gz
{ProjectName}_{SampleName}_L001_R2_001.fastq.gz
```

The sample name is extracted as everything before `_L001_R1_001`.

**Example:**
```
Project1_Sample1_L001_R1_001.fastq.gz
Project1_Sample1_L001_R2_001.fastq.gz
Project1_Sample2_L001_R1_001.fastq.gz
Project1_Sample2_L001_R2_001.fastq.gz
```

## Output Directory Structure

```
<output_dir>/
├── config.yaml                              # Pipeline config (auto-generated)
├── submit_job.sh                            # SLURM script (auto-generated)
├── 00_raw/                                  # Symlinks to input fastq files
├── 01_adapter/                              # Adapter-trimmed reads
│   ├── summary_adapter_trimming.txt
│   ├── 01_logs/
│   └── README.txt
├── 02_primer_trimmed/                       # Primer-trimmed reads
│   ├── summary_primer_trimming.txt
│   ├── 02_logs/
│   └── README.txt
├── 03_dada2/                                # Final DADA2 results
│   ├── filtered/                            # Quality-filtered fastqs
│   ├── seqtab.rds
│   ├── seqtab_nochim.rds
│   ├── sequences.fasta
│   ├── sequence_process_summary.txt
│   ├── taxID.rds
│   ├── taxID_bootstrap.rds
│   ├── <AMPLICON>_combined_sequences_taxa.txt
│   ├── <AMPLICON>_combined_sequences_taxa_bootstrap.txt
│   └── README.txt
└── logs/
    └── snakemake_TIMESTAMP.log
```

## Dry Run (Validation)

To preview the workflow without running it:
```bash
snakemake -n --snakefile Snakefile --configfile <output_dir>/config.yaml
```

To visualize the rule DAG:
```bash
snakemake --dag --configfile <output_dir>/config.yaml | dot -Tpdf > dag.pdf
```

## Parameters by Amplicon Type

### Filtering Parameters

| Amplicon | Quality=good | Quality=bad |
|----------|-------------|-----------|
| **16S-V4** | truncLen=(240,200), maxEE=(2,4) | truncLen=(240,200), maxEE=(4,6) |
| **ITS1/ITS2** | minLen=50, maxEE=(2,2) | minLen=50, maxEE=(4,4) |
| **18S-AMF** | minLen=100, maxEE=(2,4), truncQ=5, matchIDs=TRUE | minLen=100, maxEE=(4,6), truncQ=5, matchIDs=TRUE |
| **18S-V4** | minLen=100, maxEE=(2,4), truncQ=5, matchIDs=TRUE | minLen=100, maxEE=(4,6), truncQ=5, matchIDs=TRUE |

### Merging Parameters

| Amplicon | minOverlap | justConcatenate |
|----------|-----------|-----------------|
| 16S-V4 | 20 | No |
| ITS1/ITS2 | 50 | No |
| 18S-AMF | 10 | **Yes** |
| 18S-V4 | 10 | No |

#### ⚠️ 18S-V4 Read Merging

**18S-V4 datasets often have poor reverse read quality**, which can result in low merge rates. Monitor the merge rate in the QC summary:
- **Merge rate >70%**: Normal, proceed as usual
- **Merge rate 50-70%**: Acceptable, but expect some data loss
- **Merge rate <50%**: Consider using **forward reads only**

To use forward reads only, modify `workflow/scripts/run_dada2.R` line 195 to set `justConcatenate=TRUE`:
```r
merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse,
  trimOverhang=TRUE, minOverlap=10, justConcatenate=TRUE  # <-- Add this
)
```

This concatenates R1+R2 without requiring overlap, preserving all data but sacrificing merge-based quality control.

### Taxonomy Database

Stored on MSI Agate at `/projects/standard/kennedyp/shared/taxonomy/`:

| Amplicon | Database | File |
|----------|----------|------|
| 16S-V4 | SILVA | silva_nr99_v138.1_train_set.fa |
| ITS1/ITS2 | UNITE | sh_general_release_dynamic_all_19.02.2025_wKennedySynmock.fasta |
| 18S-AMF | MaarjAM | maarjam_dada2_SSU_2021.fasta |
| 18S-V4 | PR2 | pr2_version_5.1.1_SSU_dada2.fasta |

## SLURM Resource Allocation

Default allocation (suitable for most runs):
- **Time**: 12 hours
- **CPUs**: 64 cores
- **Memory**: 60 GB

Adjust in `submit_job.sh` if needed.

## Troubleshooting

### No fastq files found
Check that files match the naming pattern: `*_L001_R1_001.fastq.gz` and `*_L001_R2_001.fastq.gz`

### Conda environment fails
Ensure conda is installed and activate it before running:
```bash
source activate
```

### Permission denied on R scripts
Make sure R modules are loaded:
```bash
module load R/4.4.0-openblas-rocky8
```

### Snakemake not found
Make sure snakemake is installed in your environment or load it via module if available.

## Reference Files

The original reference scripts used to build this pipeline are in:
- `16S/` - 16S-V4 reference scripts
- `ITS1/` - ITS1 reference scripts
- `ITS2_Fun/` - ITS2 reference scripts
- `18S_AMF/` - 18S-AMF reference scripts
- `18S_MicroEuk/` - 18S-V4 reference scripts

These are for historical reference only; use the unified pipeline instead.

## Citation & Credits

Built on DADA2 (Callahan et al., 2016) and Snakemake workflows.
Designed for Kennedy Lab amplicon sequencing pipelines.

## License

Internal use only
