# Implementation Status — HTS_DADA2processing Pipeline Refactor

## Completed ✅

### Core Infrastructure
- [x] `setup.sh` — Configuration setup, ~/bin installation, R environment check
- [x] `run_dada2processing.sh` — User-facing CLI wrapper (copied to ~/bin)
- [x] `dada2processing_msiSLURM.sh` — SLURM job orchestration with routing logic
- [x] `.gitignore` — Standard git ignore patterns
- [x] `scripts/primercheck.sh` — Fixed primer detection and validation utility

### Adapter Trimming Scripts
- [x] `scripts/illumina/adapter_trimming_16S-V4.sh` — Illumina 16S-V4 adapter + primer trimming

## In Progress / To Be Completed 🔧

### Illumina Adapter Trimming Scripts (4 remaining)
- [ ] `scripts/illumina/adapter_trimming_ITS1.sh` — ITS1F/ITS2 primers
- [ ] `scripts/illumina/adapter_trimming_ITS2.sh` — 5.8SR/ITS4 primers
- [ ] `scripts/illumina/adapter_trimming_18S-V4.sh` — 616*f/1132r primers (MicroEuk)
- [ ] `scripts/illumina/adapter_trimming_18S-AMF.sh` — WANDA_1/AML2_1 primers

**Template Pattern:** Use cutadapt with `-a/-A` (adapter) and `-g/-G` (primer) flags.  
Primers for each marker gene (from UMGC defaults):
```
16S-V4:  515f (GTGYCAGCMGCCGCGGTAA) / 806r (GGACTACNVGGGTWTCTAAT)
ITS1:    ITS1F (CTTGGTCATTTAGAGGAAGTAA) / ITS2 (GCTGCGTTCTTCATCGATGC)
ITS2:    5.8SR (TCGATGAAGAACGCAGCG) / ITS4 (TCCTCCGCTTATTGATATGC)
18S-V4:  616*f (TTAAARVGYTCGTAGTYG) / 1132r (CCGTCAATTHCTTYAART)
18S-AMF: WANDA_1 (CAGCCGCGGTAATTCCAGCT) / AML2_1 (GAACCCAAACACTTTGGTTTCC)
```

### Aviti Adapter Trimming Scripts (5 required)
- [ ] `scripts/aviti/adapter_trimming_16S-V4.sh` — Adapt from `dada2processing-main/adapter_trimming_16S_aviti.sh`
- [ ] `scripts/aviti/adapter_trimming_ITS1.sh` — Adapt from `dada2processing-main/adapter_trimming_ITS.sh`
- [ ] `scripts/aviti/adapter_trimming_ITS2.sh` — Adapt from `dada2processing-main/adapter_trimming_ITS2.sh`
- [ ] `scripts/aviti/adapter_trimming_18S-V4.sh` — Adapt from `dada2processing-main/adapter_trimming_18S_MicroEuk.sh`
- [ ] `scripts/aviti/adapter_trimming_18S-AMF.sh` — Adapt from `dada2processing-main/adapter_trimming_18S_AMF.sh`

**Note:** These require `scripts/aviti/aviti_adapters.fa` and `scripts/aviti/aviti_ITSprimer_as_adapter.fa`  
→ **Pinned:** Request from Trevor Gould. For now, scripts should gracefully handle missing files.

### Illumina DADA2 R Scripts (4 required)
- [ ] `scripts/illumina/run_16S-V4_dada2.R` — Adapt from `dada2processing-main/run_16S_dada2.R`
  - Accept `args[2]` = output_dir
  - Taxonomy: SILVA v138 from Kennedy lab shared folder
  - Quality logic: keep as-is (good/bad modes)

- [ ] `scripts/illumina/run_ITS_dada2.R` — Adapt from `dada2processing-main/run_ITS_dada2.R`
  - Single script handles both ITS1 and ITS2 (primer trimming is separate)
  - Accept `args[2]` = output_dir
  - Taxonomy: UNITE from Kennedy lab shared folder
  - Quality logic: keep as-is

- [ ] `scripts/illumina/run_18S-V4_dada2.R` — **NEW**: Create based on `run_18S_dada2.R`
  - Target: Protists (MicroEuk) with 18S-V4 primers
  - Taxonomy: PR2 from Kennedy lab shared folder
  - Quality logic: keep same as run_18S_dada2.R but only use MaarjAM database (drop EUK reference)

- [ ] `scripts/illumina/run_18S-AMF_dada2.R` — Adapt from `dada2processing-main/run_18S_dada2.R`
  - Target: AMF (Arbuscular Mycorrhizal Fungi)
  - Taxonomy: MaarjAM **only** (drop DADA2_EUK reference)
  - Quality logic: keep as-is
  - Accept `args[2]` = output_dir

**Common Changes for all Illumina R scripts:**
- Add `output_dir <- args[2]` at line ~18
- Replace all `"../dada2output/"` with `output_dir`
- Update taxonomy paths to `/panfs/jay/groups/4/kennedyp/shared/taxonomy/`

### Aviti DADA2 R Scripts (4 required)
- [ ] `scripts/aviti/run_16S-V4_dada2.R` — Adapt from `dada2processing-main/run_16S_dada2_aviti.R`
  - **FIX:** Re-enable quality branching (currently all commented out)
  - Add `args[2]` = output_dir
  - Taxonomy: SILVA from Kennedy lab shared folder
  - Update paths

- [ ] `scripts/aviti/run_ITS_dada2.R` — Adapt from `dada2processing-main/run_ITS_dada2_aviti.R`
  - **FIX:** Update UNITE database from April 2024 → February 2025 release
  - **FIX:** Re-enable quality branching (currently commented out)
  - Add `args[2]` = output_dir
  - Taxonomy: UNITE from Kennedy lab shared folder

- [ ] `scripts/aviti/run_18S-V4_dada2.R` — **NEW**: Create based on Aviti 18S template
  - Taxonomy: PR2 from Kennedy lab shared folder
  - Quality logic & error modeling: same as run_18S_dada2_aviti.R
  - Add `args[2]` = output_dir

- [ ] `scripts/aviti/run_18S-AMF_dada2.R` — Adapt from `dada2processing-main/run_18S_dada2_aviti.R`
  - Keep: MaarjAM taxonomy only
  - Keep: binned error correction if present
  - Add `args[2]` = output_dir
  - Update paths & database references

**Common Changes for all Aviti R scripts:**
- Add `output_dir <- args[2]` at line ~18
- Replace all `"../dada2output/"` with `output_dir`
- Update taxonomy paths to `/panfs/jay/groups/4/kennedyp/shared/taxonomy/`

### Documentation
- [ ] `README.md` — Installation, usage, examples, troubleshooting

## Architecture Notes

### Directory Flow
1. **User input:** Raw fastq files in `<input_dir>`
2. **After adapter trimming:** Trimmed paired reads in `01_adapter/` (or intermediate step)
3. **After primer trimming:** Filtered reads in `02_filtered/`
4. **After DADA2:** ASV table + taxonomy in `<output_dir>/`
5. **Logs:** Summary logs in `<output_dir>/logs/` or `<output_dir>/`

### Taxonomy Database Paths (Kennedy Lab Shared Folder)
Expected path structure: `/panfs/jay/groups/4/kennedyp/shared/taxonomy/`
- SILVA: `silva_*` files
- UNITE: `UNITE_*` or similar
- PR2: `PR2_*` files
- MaarjAM: `maarjam_*` files

These paths should match what HTS_ASV2OTU uses.

## Testing Checklist

Once scripts are created:
1. [ ] Run `bash setup.sh` (verify config.sh is created, ~/bin/run_dada2processing is installed)
2. [ ] Verify `run_dada2processing --help` works from any directory
3. [ ] Test with small sample dataset: `run_dada2processing /path/to/test_fastq /path/to/output 16S-V4 illumina`
4. [ ] Check `tail -f <output_dir>/pipeline_*.out` for job progress
5. [ ] Verify output files are created in `<output_dir>/`
6. [ ] Test different marker genes and platforms

## Files in dada2processing-main/ for Reference

Source scripts to adapt (in order of priority):
```
dada2processing-main/
├── run_16S_dada2.R                (→ scripts/illumina/run_16S-V4_dada2.R)
├── run_16S_dada2_aviti.R          (→ scripts/aviti/run_16S-V4_dada2.R)
├── run_ITS_dada2.R                (→ scripts/illumina/run_ITS_dada2.R)
├── run_ITS_dada2_aviti.R          (→ scripts/aviti/run_ITS_dada2.R)
├── run_18S_dada2.R                (→ scripts/illumina/run_18S-V4/AMF_dada2.R)
├── run_18S_dada2_aviti.R          (→ scripts/aviti/run_18S-V4/AMF_dada2.R)
├── adapter_removal.sh             (→ scripts/illumina/adapter_trimming_16S-V4.sh) ✅
├── adapter_trimming_16S_aviti.sh  (→ scripts/aviti/adapter_trimming_16S-V4.sh)
├── adapter_trimming_ITS.sh        (→ scripts/aviti/adapter_trimming_ITS1.sh)
├── adapter_trimming_ITS2.sh       (→ scripts/aviti/adapter_trimming_ITS2.sh)
├── adapter_trimming_18S_MicroEuk.sh (→ scripts/aviti/adapter_trimming_18S-V4.sh)
└── adapter_trimming_18S_AMF.sh    (→ scripts/aviti/adapter_trimming_18S-AMF.sh)
```
