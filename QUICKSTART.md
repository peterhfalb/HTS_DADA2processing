# DADA2 Pipeline Quick Start

## Setup (One-time)

```bash
# Navigate to the pipeline directory
cd /path/to/DADA2_Snakemake_groundUp

# Run the installer
bash install.sh

# Reload your shell
source ~/.bashrc

# Verify installation
run_dada2processing --help
```

## Running the Pipeline

### Example: 16S-V4 Analysis

```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir /path/to/raw/fastq/files \
  --output-dir /path/to/output/directory \
  --email your.email@umn.edu
```

### Available Amplicon Types
- `16S-V4` - Bacteria/Archaea
- `ITS1` - Fungi (ITS1 region)
- `ITS2` - Fungi (ITS2 region)
- `18S-AMF` - Arbuscular mycorrhizal fungi
- `18S-V4` - Microbial eukaryotes

### Input Files
Files in `--fastq-dir` must be named:
```
{ProjectName}_{SampleName}_L001_R1_001.fastq.gz
{ProjectName}_{SampleName}_L001_R2_001.fastq.gz
```

Example:
```
Project1_Sample1_L001_R1_001.fastq.gz
Project1_Sample1_L001_R2_001.fastq.gz
Project1_Sample2_L001_R1_001.fastq.gz
Project1_Sample2_L001_R2_001.fastq.gz
```

## Job Submission

After running `run_dada2processing`, a SLURM script is auto-generated:

```bash
# Submit the job to Agate
sbatch <output_dir>/submit_job.sh

# Monitor progress
tail -f <output_dir>/logs/snakemake_*.log
```

## Advanced Options

### Quality Setting
Use stricter filtering for high-quality reads, or more lenient for low-quality:
```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir ./fastqs \
  --output-dir ./output \
  --email you@umn.edu \
  --quality bad
```

Options: `good` (default), `bad`

### Parallel Jobs
Adjust max concurrent Snakemake jobs (default: 16):
```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir ./fastqs \
  --output-dir ./output \
  --email you@umn.edu \
  --jobs 8
```

## Output Files

Key results in `<output_dir>/03_dada2/`:

| File | Description |
|------|-------------|
| `{AMPLICON}_combined_sequences_taxa.txt` | **Main output**: ASV table with taxonomy |
| `{AMPLICON}_combined_sequences_taxa_bootstrap.txt` | With bootstrap confidence scores |
| `sequences.fasta` | Fasta file of all unique sequences |
| `sequence_process_summary.txt` | Per-sample processing statistics |
| `seqtab_nochim.rds` | R object: final sequence abundance table |

## Troubleshooting

### "No fastq files found"
- Check file naming pattern (must have `_L001_R1_001.fastq.gz` and `_L001_R2_001.fastq.gz`)
- Verify directory is correct

### SLURM job fails to start
- Check `<output_dir>/logs/snakemake_*.log` for errors
- Verify you have write permissions in `--output-dir`
- Check that Agate has the required modules: `module avail R/4.4` and `module avail cutadapt`

### Conda environment issues
- Ensure you're on Agate (not a login node)
- Try: `module load conda` or ensure conda is in your PATH

## Tips

### Dry-run (no execution)
Preview the workflow without running:
```bash
snakemake -n --snakefile Snakefile --configfile <output_dir>/config.yaml
```

### Visualize workflow
Create a DAG visualization:
```bash
snakemake --dag --configfile <output_dir>/config.yaml | dot -Tpdf > dag.pdf
```

### Re-run failed jobs
To rerun from where it failed (Snakemake auto-detects incomplete outputs):
```bash
sbatch <output_dir>/submit_job.sh
```

## Questions?

Refer to the full README.md for detailed documentation.
