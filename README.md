# HTS_DADA2processing — Kennedy Lab Microbiome Pipeline

A containerized, production-ready DADA2 processing pipeline for microbiome amplicon sequencing data. Processes raw FASTQ files through adapter/primer trimming and DADA2 denoising to generate amplicon sequence variants (ASVs) with taxonomy assignments.

**Maintained by:** Kennedy Lab, University of Minnesota  
**Original pipeline by:** Trevor Gould  
**Refactored for Kennedy Lab:** March 2026

---

## Features

- ✅ **Single command execution:** `run_dada2processing` — easy-to-use CLI
- ✅ **SLURM-integrated:** Runs on MSI AGATE supercomputer
- ✅ **Multi-marker support:** 16S-V4, ITS1, ITS2, 18S-V4 (protists), 18S-AMF
- ✅ **Multi-platform support:** Illumina (MiSeq/NovaSeq) and Aviti (Element Biosciences)
- ✅ **Automated primer validation:** Detects and warns on primer mismatches
- ✅ **Flexible quality control:** "good" (default) or "bad" quality filtering modes
- ✅ **Standardized output:** Consistent taxonomy tables across all marker genes

---

## Installation

### 1. Clone the Repository

```bash
ssh -Y <umn_username>@agate.msi.umn.edu
mkdir -p ~/pipelines
cd ~/pipelines
git clone https://github.com/peterhfalb/HTS_DADA2processing.git
cd HTS_DADA2processing
```

### 2. Run Setup

```bash
bash setup.sh
```

When prompted, enter your UMN email address for SLURM notifications.

### 3. Update Your PATH

The setup script installs the `run_dada2processing` command to `~/bin`. Ensure this is in your PATH:

```bash
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### 4. Verify Installation

```bash
run_dada2processing --help
```

### 5. (Future) Request Aviti Adapter Files from Trevor Gould

The pipeline requires two FASTA files for Aviti adapter trimming. Contact Trevor Gould to obtain:
- `aviti_adapters.fa`
- `aviti_ITSprimer_as_adapter.fa`

Place them in `scripts/aviti/` directory:

```bash
cp aviti_adapters.fa ~/pipelines/HTS_DADA2processing/scripts/aviti/
cp aviti_ITSprimer_as_adapter.fa ~/pipelines/HTS_DADA2processing/scripts/aviti/
```

---

## Usage

### Command Syntax

```bash
run_dada2processing <input_fastq_dir> <output_dir> <marker_gene> <platform> [--quality bad]
```

### Arguments

| Argument | Required | Options | Description |
|---|---|---|---|
| `input_fastq_dir` | yes | path | Directory containing raw `_R1_001.fastq.gz` and `_R2_001.fastq.gz` files |
| `output_dir` | yes | path | Output directory (will be created if needed) |
| `marker_gene` | yes | `16S-V4`, `ITS1`, `ITS2`, `18S-V4`, `18S-AMF` | Sequencing target |
| `platform` | yes | `illumina`, `aviti` | Sequencing platform |
| `--quality` | optional | `good` (default), `bad` | Filtering stringency |

### Examples

```bash
# 16S-V4 Illumina (good quality — default)
run_dada2processing /path/to/raw/fastq /path/to/output 16S-V4 illumina

# ITS2 Aviti (bad quality)
run_dada2processing /path/to/raw/fastq /path/to/output ITS2 aviti --quality bad

# 18S-V4 protists (Illumina)
run_dada2processing /path/to/raw/fastq /path/to/output 18S-V4 illumina

# 18S-AMF (Aviti)
run_dada2processing /path/to/raw/fastq /path/to/output 18S-AMF aviti
```

### Monitoring Progress

Once submitted, monitor the SLURM job:

```bash
# View real-time output
tail -f <output_dir>/pipeline_*.out

# Check job status
squeue -u $USER
```

---

## Quality Modes Explained

| Mode | Forward maxEE | Reverse maxEE | Forward truncLen | Reverse truncLen | Use When |
|---|---|---|---|---|---|
| **good** | 2 | 2 | 240 (16S) | 160 (16S) | High-quality sequencing, strict filtering |
| **bad** | 4 | 6 | 240 (16S) | 160 (16S) | Lower-quality sequencing, relaxed filtering |

### Choosing Quality Mode

- **Default to "good"** for typical Illumina/Aviti runs
- Use **"bad"** if:
  - Quality scores are visibly low (FastQC shows many low-quality bases)
  - You're losing too many reads with "good" mode
  - Your sequencing came from a challenging sample (old DNA, environmental samples, etc.)

---

## Marker Genes & Databases

| Marker Gene | Primers | Primers (F/R) | Database | Notes |
|---|---|---|---|---|
| **16S-V4** | 515f/806r | GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT | SILVA v138 | Bacteria & archaea |
| **ITS1** | ITS1F/ITS2 | CTTGGTCATTTAGAGGAAGTAA / GCTGCGTTCTTCATCGATGC | UNITE v9 | Fungi |
| **ITS2** | 5.8SR/ITS4 | TCGATGAAGAACGCAGCG / TCCTCCGCTTATTGATATGC | UNITE v9 | Fungi |
| **18S-V4** | 616f/1132r | TTAAARVGYTCGTAGTYG / CCGTCAATTHCTTYAART | PR2 v4.14 | Protists (microeukaryotes) |
| **18S-AMF** | WANDA_1/AML2_1 | CAGCCGCGGTAATTCCAGCT / GAACCCAAACACTTTGGTTTCC | MaarjAM v2.40 | Arbuscular mycorrhizal fungi |

All databases are stored in Kennedy Lab shared taxonomy folder: `/projects/standard/kennedyp/shared/taxonomy/`

---

## Output Files

After successful completion, the output directory will contain:

### ASV Table (Primary Output)
- `<marker>_combined_sequences_taxa_<database>.txt` — OTU abundance table with taxonomy
- `<marker>_combined_sequences_taxa_<database>_boot.txt` — Same with bootstrap confidence

### Supporting Files
- `sequences.fasta` — ASV sequences in FASTA format
- `seqtab.rds` / `seqtab_nochim.rds` — R objects for downstream analysis
- `sequence_process_summary.txt` — Statistics on reads lost at each step
- `taxID<database>.rds` — Taxonomy assignments (R format)
- `taxID<database>_bootstrap.rds` — Bootstrap confidence values (R format)
- `pipeline_*.out` / `pipeline_*.err` — SLURM job logs

### Directory Structure After Completion

```
<output_dir>/
├── 16S-V4_combined_sequences_taxa_silva.txt       (primary output)
├── 16S-V4_combined_sequences_taxa_silva_boot.txt
├── sequences.fasta
├── sequence_process_summary.txt
├── seqtab.rds
├── seqtab_nochim.rds
├── taxIDsilva.rds
├── taxIDsilva_bootstrap.rds
├── pipeline_JOBID.out                             (SLURM log)
├── logs/                                          (adapter/primer trimming logs)
├── 01_adapter/                                    (intermediate: after adapter removal)
└── 02_filtered/                                   (intermediate: after primer trimming)
```

---

## Pipeline Architecture

```
Raw FASTQ files
        ↓
   Primer Detection (primercheck.sh)
        ↓
   Adapter/Primer Trimming (platform-specific)
        ↓
   DADA2 Processing (R script)
        ├─ Quality filtering
        ├─ Dereplication
        ├─ Error learning
        ├─ ASV inference
        ├─ Pair merging
        └─ Chimera removal
        ↓
   Taxonomy Assignment (SILVA/UNITE/PR2/MaarjAM)
        ↓
   Output: ASV table with taxonomy
```

### File Naming Conventions

- **Raw reads:** `SAMPLENAME_R1_001.fastq.gz`, `SAMPLENAME_R2_001.fastq.gz`
- **Intermediate (adapter-trimmed):** `01_adapter/SAMPLENAME_R1_001.fastq.gz`
- **Final (primer-trimmed):** `02_filtered/SAMPLENAME_R1_001.fastq.gz`

---

## Troubleshooting

### Common Issues

#### 1. "config.sh not found"
**Solution:** Run `bash setup.sh` in the HTS_DADA2processing directory first.

#### 2. "run_dada2processing command not found"
**Solution:** Ensure `~/bin` is in your PATH:
```bash
echo $PATH | grep "$HOME/bin"
echo 'export PATH=$HOME/bin:$PATH' >> ~/.bashrc && source ~/.bashrc
```

#### 3. "Aviti adapter FASTA files not found"
**Solution:** Request `aviti_adapters.fa` and `aviti_ITSprimer_as_adapter.fa` from Trevor Gould, then place in `scripts/aviti/`.

#### 4. "primer check shows 0 reads" 
**Solution:** This usually means primers have already been removed (expected). The pipeline will continue normally.

#### 5. "No fastq.gz files found in input directory"
**Solution:** Verify files exist and have the correct naming pattern:
```bash
ls <input_dir>/*_R1_*.fastq.gz
ls <input_dir>/*_R2_*.fastq.gz
```

#### 6. Job times out or runs out of memory
**Solution:** The job defaults to 12 hours, 64 CPUs, 60 GB RAM. For larger datasets, contact the Kennedy Lab to adjust SLURM parameters.

---

## Citation

If you use this pipeline, please cite:

1. **DADA2:** Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicons. *Nat Methods*. DOI: 10.1038/nmeth.3869

2. **Taxonomy databases:**
   - SILVA: Quast et al. (2013) — https://www.arb-silva.de/
   - UNITE: Nilsson et al. (2019) — https://unite.ut.ee/
   - PR2: Guillou et al. (2013) — https://pr2-database.org/
   - MaarjAM: Opik et al. (2010) — https://maarjam.botany.ut.ee/

---

## Support & Questions

For issues or questions about the pipeline:
1. Check the troubleshooting section above
2. Review SLURM output logs in the output directory
3. Contact the Kennedy Lab for configuration or run-specific questions

---

## Version History

- **March 2026:** Initial refactor for Kennedy Lab (Illumina + Aviti support, multi-marker)
- **Original:** Trevor Gould's DADA2 processing scripts (2020-2025)

---

## License

This pipeline builds on Trevor Gould's original work. Please see the LICENSE file for details.

