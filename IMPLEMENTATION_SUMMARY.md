# DADA2 Snakemake Pipeline - Implementation Summary

## Project Completion

The DADA2 Snakemake pipeline has been successfully implemented with support for all five amplicon types and flexible quality/parameter controls.

## Files Created

### Core Pipeline
- **`Snakefile`** - Main workflow orchestration
- **`bin/run_dada2processing`** - CLI entry point (must be installed to PATH)
- **`install.sh`** - Setup script to add pipeline to PATH

### Workflow Rules
- **`workflow/rules/cutadapt.smk`** - Adapter and primer trimming rules
- **`workflow/rules/dada2.smk`** - DADA2 denoising and taxonomy rules

### Scripts
- **`workflow/scripts/run_dada2.R`** - Unified parameterized R script for all amplicon types

### Configuration
- **`envs/dada2.yaml`** - Conda environment specification

### Documentation
- **`README.md`** - Comprehensive user guide
- **`IMPLEMENTATION_SUMMARY.md`** - This file

## Key Features

### 1. CLI-Based Invocation
```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir /path/to/raw \
  --output-dir /path/to/output \
  --email user@umn.edu \
  [--quality good|bad] \
  [--jobs 16]
```

- Single command interface
- Automatic sample discovery from file names
- Config file auto-generation
- SLURM script auto-generation

### 2. Amplicon Type Support
- **16S-V4**: Bacteria/Archaea (SILVA database)
- **ITS1**: Fungi ITS1 (UNITE database)
- **ITS2**: Fungi ITS2 (UNITE database)
- **18S-AMF**: Arbuscular mycorrhizae (MaarjAM database)
- **18S-V4**: Microbial eukaryotes (PR2 database)

### 3. Quality-Based Parameterization
Each amplicon type has two parameter sets:
- **`--quality good`** (default): Stricter thresholds (lower maxEE values)
- **`--quality bad`**: More permissive thresholds (higher maxEE values)

### 4. Processing Pipeline
Three-stage architecture:
1. **Adapter trimming** (cutadapt, 8 cores per sample)
2. **Primer trimming** (cutadapt, 4 cores per sample)
3. **DADA2 processing** (filtering, denoising, merging, taxonomy, single job)

### 5. Output Organization
Structured output with READMEs at each stage:
- `00_raw/` - Input symlinks
- `01_adapter/` - Adapter-trimmed reads + logs
- `02_primer_trimmed/` - Primer-trimmed reads + logs
- `03_dada2/` - Final results (seqtab, fasta, taxonomy)

## Technical Highlights

### Snakemake Integration
- Declarative workflow definition
- Automatic parallelization of per-sample rules
- Conda environment management
- Dry-run support for validation

### Parameter Configuration
Primer sequences and taxonomy databases are defined in the Snakefile as Python dicts and automatically selected based on `--amplicon`:
```python
PRIMERS = {
    "16S-V4": {"fwd": "...", "rev": "..."},
    ...
}
TAXONOMY_DB = {
    "16S-V4": "/projects/standard/...",
    ...
}
```

### Unified R Script
Single `run_dada2.R` script handles all amplicon types with conditional logic:
- Different filterAndTrim parameters per amplicon + quality
- Different mergePairs parameters per amplicon
- Length filtering for 18S amplicons only
- Amplicon-specific taxonomyAssignment parameters

### SLURM Integration
- Auto-generated SLURM script template
- Default resources: 12h, 64 cores, 60GB
- Email notifications
- Easy customization before submission

## Installation & Usage

### 1. Install
```bash
cd /path/to/DADA2_Snakemake_groundUp
bash install.sh
source ~/.bashrc
```

### 2. Run
```bash
run_dada2processing \
  --amplicon 16S-V4 \
  --fastq-dir ./fastqs \
  --output-dir ./output \
  --email you@umn.edu
```

### 3. Submit
```bash
sbatch <output_dir>/submit_job.sh
```

### 4. Monitor
```bash
tail -f <output_dir>/logs/snakemake_*.log
```

## Design Decisions

### Single R Script vs. Multiple Scripts
**Choice**: Single unified R script with conditional logic
**Rationale**: Reduces code duplication, easier to maintain parameter differences, single entry point

### Config File Generation vs. Template
**Choice**: Generate fresh config.yaml per run
**Rationale**: Simpler user experience, auto-detection of samples, no template maintenance needed

### Primer/DB Lookup in Snakefile vs. External Config
**Choice**: Python dicts in Snakefile
**Rationale**: Primer sets and DB paths are static and tightly coupled to amplicon types; simplifies the config file and reduces chance of mismatch

### Output Directory Structure
**Choice**: Numbered prefixes (01_, 02_, 03_) with separate subdirectories for logs
**Rationale**: Clear progression through pipeline stages, easy to navigate, logs isolated from data

## Testing Recommendations

### Dry Run
```bash
snakemake -n --snakefile Snakefile --configfile <output_dir>/config.yaml
```

### Visualize DAG
```bash
snakemake --dag --configfile <output_dir>/config.yaml | dot -Tpdf > dag.pdf
```

### End-to-End Test
1. Use existing test fastq files from one of the reference folders
2. Run pipeline with `--quality good` first
3. Verify all output files are present
4. Check that summary statistics look reasonable

## Known Limitations & Future Improvements

### Current Limitations
- Taxonomy database paths are MSI-specific; must be updated for other clusters
- R/4.4.0-openblas-rocky8 module load is MSI-specific

### Possible Future Enhancements
- Externalize MSI paths to a config file for portability
- Add rule for ASV table filtering/decontamination
- Add OTU clustering option as alternative to ASV
- Add rarefaction/diversity analysis rules
- Generate HTML report of processing statistics

## References

- DADA2: Callahan et al. (2016), Nature Methods
- Snakemake: Mölder et al. (2021), F1000Research
- Taxonomic databases: SILVA, UNITE, PR2, MaarjAM
