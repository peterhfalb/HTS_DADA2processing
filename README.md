# SLURM pipeline to convert demultiplexed fastq files to denoised ASV tables and denoised OTU tables

### Pipeline Overview and Authors

This pipeline script was written specifically for use by the Kennedy Lab at the University of Minnesota. It takes demultiplexed fastq files coming from the UMGC sequencing center, trims adapters and primers, then runs denoising using the DADA2 pipeline. Then, ASV tables are generated using the DADA2 taxonomy assigner. The pipeline then does further denoising using an OTU-based pipeline to further reduce sequencing noise and errors coming off the sequencing platform, and will likely be most helpful for cleaning non-host associated bulk communities with extremely high sequence diversity. OTUs are clustered at 97% sequence identity, which very approximately represents the level of species. The pipeline is specifically designed for the Minnesota Agate computing cluster, and requires user access to the Kennedy lab shared folders for use. If you would like to run this on a different computing cluster, or don't have access to the Kennedy Lab folders, please contact me (Peter F; falb0011@umn.edu). 

Adapter trimming, and DADA2 denoising settings were developed by Trevor Gould and adapted for use with this pipeline. Credit for the original development of the ASV to OTU pipeline goes to Eivind Ronold and Luis Morgado at the University of Oslo. I (PF) have updated and adapted all the scripts to work on University of Minnesota systems, as well as implemented various bug fixes and speed/code improvements. Individual scripts were consolidated into one pipeline script which compiles and installs all necessary packages for running the pipeline. The pipeline was also adapted to handle all types of sequencing data currently used in the Kennedy Lab (16S-V4, ITS1, ITS2, 18S-V4, 18S-AMF). 

## Installation instructions

Login in to MSI via terminal on your mac (usage might be different on Windows, see https://msi.umn.edu/connecting/connecting-to-hpc-resources):

```bash
# login to MSI via SSH secure shell
ssh -Y yourUMNusername@agate.msi.umn.edu # change 'yourUMNusername' to your UMN username
```

Next, move into a personal directory where you would like the pipeline files to be installed, or create a new directory to house your pipelines. This is a permanent, single location you will use anytime you run the pipeline, and is personal to you (it won't work for other users).

```bash
# [optional]: make a directory for pipelines/packages (if desired/needed)
mkdir pipelines # mkdir is "make directory"

# move into the directory where you want the pipeline installed
cd pipelines/ # cd is "change directory"
```

Once in your preferred directory, clone the HTS_ASV2OTU github repository:
```bash
git clone https://github.com/peterhfalb/HTS_DADA2processing.git
```

To set up and install the pipeline, open the cloned repository and run the setup/installation script:
```bash
# move into the cloned repository
cd HTS_DADA2processing/ 

# run this line of code to install and compile the pipeline package
bash install.sh
source ~/.bashrc
```
**NOTE: It may take 10-15 minutes to install the necessary packages**

This install script will create a command you can run from anywhere called `run_dada2processing`

After proper installation, the DADA2 processing pipeline should now be ready to go.

## Updating the Pipeline:

As improvements and bug fixes are made to the pipeline, you can easily update your local installation using git. From your pipeline directory, simply run:

```bash
cd ~/your/path/to/HTS_DADA2processing
git pull
```

After a git pull it is usually best practice to reinstall the virtual environment, just in case any package dependencies have changed:

```bash
# run the install script
bash install.sh
# when prompted, say 'y' to delete and reinstall the virtual environment
```
This will download the latest version and update all pipeline scripts.

**If git pull gives you errors:**

If you encounter git-related errors when trying to pull updates (e.g., merge conflicts, local changes preventing pull), the simplest solution is to delete the pipeline folder and re-clone:

```bash
# Navigate to the parent directory
cd ~/your/path/to/
rm -rf HTS_DADA2processing

# Re-clone the repository
git clone https://github.com/peterhfalb/HTS_DADA2processing.git

# Navigate into the cloned repository
cd HTS_DADA2processing

# Re-run setup to ensure everything is configured
bash install.sh
```

Your `run_dada2processing` command and conda environment will be reconfigured and ready to use.

## Usage Instructions:

This script takes one main input -- a folder containing your demultiplexed fastq files. You will also need to specify an output folder on the cluster, which should be different than the fastq folder. 

**Note:** *You can upload/delete files and folders using an interactive GUI if that is easier. To access that, go to msi.umn.edu, then click "Open OnDemand", then click "Files" and go to either your Home Directory, or the shared kennedyp folder. Now, you can click through your folders, and edit as desired. No computational tasks can be run through here, however.*

**Step 1: Login to the cluster, and prepare folder**
```bash
# We need two things: the path to your input fastq files, and the path the the output folder
# FIRST: locate where your fastq files are and take note of the folder path
# SECOND: locate or create a folder for your output files (see below):

# login to the cluster (fill in your UMN username)
ssh -Y yourMSIusername@agate.msi.umn.edu 

# navigate to where you want your project to live
cd yourPreferredDirectory/ 

# create a directory for your project
mkdir yourProjectName_ASVtoOTU 

# navigate into the new directory
cd yourProjectName_ASVtoOTU 

# print and copy the path to the output file directory
pwd # "print working directory" - copy what this prints

# leave the ssh window
exit 
```

**Step 2: Run the pipeline**

You can view full usage information at any time by running:
```bash
run_dada2processing --help
```

To run the script, ssh login to the cluster and run `run_dada2processing` from anywhere — no need to navigate to the pipeline directory. You will have to specify at least 6 things:

*Required:*
1. Path to input fastq file folder (`--fastq-dir`)
2. Path to output folder (`--output-dir`)
3. Project Name, this will determine how your output files are named (`--project-name`)
4. Amplicon / Barcoding Region (`--amplicon`)
5. Sequencing platform (`--platform`)
6. Email for SLURM submissions (`--email`)

*Optional:*
1. Sequence quality (`--quality`)
        Right now, the default for sequence quality is `good`, by setting to `bad`, some of the filtering parameters in DADA2 become more stringent. Reach out to Trevor Gould or the UMGC for further information as to whether you sequencing run is *good* or *bad* quality.
2. Taxonomy database - if using something other than default (`-taxonomy-database`)
        Manually specify a reference database for taxonomy assignment other than default. Mostly relevant for 18S-V4 and 18S-AMF (and possibly AMF). Defaults is usually the best.
3. Forward/Reverse Primer Sequence (`--fwd-primer/--rev-primer`)
        Use to manually input a primer sequence if either a. the pipeline detects a different primer than default or b. you used a non-canonical primer not included on the Kennedy Lab primer list.
4. Whether to skip OTUing (`--skip-otu`)
        OTUing is generally recommended under the following scenarios: bulk environmental non-host-associated soil samples, Aviti sequencing, non-bacterial primer sets (ITS, 18S). OTUing is however generally applicable across all datasets and will serve to further reduce possible noise in the dataset.
5. Whether to run ITSx (`--run-ITSx`)
        *ITSx is a program which removes the highly conserved regions flanking the ITS variable region. While this can improve clustering by preventing inflated similarity from conserved regions, ITSx can also remove or truncate sequences with non-standard ITS structures, including synthetic mock community members. **By default, ITSx is disabled to preserve full-length ITS sequences and avoid issues with mock communities.** If your dataset does not include synthetic controls and you want the benefits of ITS region extraction, use the `--run-itsx` flag to enable it.*    
6. Manually specify cluster percentage level (`--otu-cluster-id`)
        Clustering is by default set to 97% which is the most commonly accepted species-level threshold across groups. However, if you want to be more or less conservative, you can increase or decrease the similarity threshold. The bacterial SSU gene (16S-V4) is generally more conserved, so there may an argument for raising the similarity threshold for that group.
7. Manually specify MUMU blast ID threshold (`--mumu-blast-id`)
        The MUMU algorithm requires a minimum similarity threshold for sequences to be merged as mother-daughter. This default is 94% for bacterial datasets, and 84% for all others. Increase the default value if you want MUMU merges to be more stringent/conservative.
8. Manually specify MUMU minimum ratio (`--mumu-ratio`)
        The MUMU algorithm requires a minimum average ratio between parent-daughter abundances for them to be accepted as a merge. By default it is 100 for 16S and 1 for all other datasets. Increase the minimum ratio if you want MUMU merges to be more stringent/conservative.
9. Keep only forward reads (`--fwd-reads-only`)
        You can discard you reverse reads if merging is failing/discarding many of your reads due to lack of overlap or poor R2 quality. This will usually not be relevant but may be desired for 18S-V4 datasets.

Expect a runtime of 10 - 120 minutes. It should not go much longer than that, but the SLURM script requests 12 hours of time on the cluster just in case.

**Run the pipeline from anywhere**
```bash
# login to the cluster (fill in your UMN username)
ssh -Y yourMSIusername@agate.msi.umn.edu 

# submit the slurm job with the following command (run from anywhere):
run_dada2processing [OPTIONS]

# Required arguments:
#  --amplicon AMPLICON       Amplicon type: 16S-V4, ITS1, ITS2, 18S-AMF, 18S-V4
#  --fastq-dir DIR           Directory containing raw demultiplexed fastq files
#  --output-dir DIR          Output directory (will be created if needed)
#  --email EMAIL             Email for SLURM notifications
#  --platform PLATFORM       Sequencing platform: illumina or aviti
#  --project-name NAME       Project name for output file naming

# Optional arguments:
#  --quality QUALITY         Sequence quality level: good (default) or bad
#  --taxonomy-database DB    Override taxonomy database: SILVA, UNITE, PR2, MaarjAM, EukaryomeSSU, EukaryomeITS
#  --fwd-primer SEQUENCE     Override forward primer for this amplicon (DNA string)
#  --rev-primer SEQUENCE     Override reverse primer for this amplicon (DNA string)
#  --skip-otu                Skip OTU clustering pipeline (default: run OTU pipeline)
#  --run-itsx                Run ITSx extraction (ITS1/ITS2 only; default: off)
#  --otu-cluster-id ID       VSEARCH clustering identity (default: 0.97)
#  --mumu-blast-id PCT       BLAST identity % for mumu (default: 94 for 16S-V4, 84 for others)
#  --mumu-ratio RATIO        mumu minimum_ratio (default: 100 for 16S-V4, 1 for others)
#  --jobs JOBS               Max parallel Snakemake jobs (default: 16)
#  --help                    Show this help message

# Examples:
#  run_dada2processing --amplicon 16S-V4 --fastq-dir /path/to/raw_fastqs --output-dir /path/to/FAB2_16S-V4_output --email name0023@umn.edu --platform illumina --project-name MyProject
#  run_dada2processing --amplicon ITS2 --fastq-dir /path/to/raw_fastqs --output-dir /path/to/output --email me@umn.edu --platform aviti --project-name ITS_Study --quality bad
#  run_dada2processing --amplicon 18S-V4 --fastq-dir /path/to/raw_fastqs --output-dir /path/to/output --email me@umn.edu --platform illumina --project-name MyProject --taxonomy-database EukaryomeSSU
#  run_dada2processing --amplicon ITS2 --fastq-dir /path/to/raw_fastqs --output-dir /path/to/output --email me@umn.edu --platform illumina --project-name MyProject --fwd-primer GTGARTCATCGAATCTTTG
```

After you run the command, the pipeline will do some pre-submission checks, generate a config file with all your presets, and then validate that the primers you chose (default unless otherwise specified) are actually present in your forward and reverse reads. If it detects them in high abundance it will tell you all is good! If not, it will show you what was detected, and either automatically select another primer (if it found it at high enough abundance), or recommend that you review your files and primer choices and manually specify a different primer.

If the primer check looks good and all the other settings look correct, then you can proceed to submit the job. The end of the script will provide you with all the information, and provide an *"sbatch"* command for you to copy-paste to submit the job. It will look something like this:

```bash
SUMMARY & RECOMMENDATIONS
==================================================
✓ Both expected primers detected at >80% frequency. Good!

Do NOT proceed with SLURM submission until primers are confirmed!

[Primer detection results saved]
Config written to: /users/9/user0011/ProjectName_AMF_output/config.yaml

==========================================
DADA2 Pipeline Ready
==========================================
--SETTINGS AND SPECIFICATIONS--
Amplicon:         18S-AMF
Quality:          good
Platform:         aviti
Samples:          3
Project name:     ProjectName
Taxonomy DB:      MaarjAM (default)

--PRIMERS--
Fwd primer:       CAGCCGCGGTAATTCCAGCT
Rev primer:       GAACCCAAACACTTTGGTTTCC

--ASV2OTU SETTINGS--
OTU clustering:   0.97 identity
ITSx extraction:  disabled
mumu BLAST ID:    84%
mumu ratio:       1

--FILES AND DIRECTORIES--
Input directory:  /users/9/user0011/ProjectName_AMF_output/fastq
Output directory: /users/9/user0011/ProjectName_AMF_output
Config file:      /users/9/user0011/ProjectName_AMF_output/config.yaml
SLURM script:     /users/9/user0011/ProjectName_AMF_output/submit_job.sh
Log file:         /users/9/user0011/ProjectName_AMF_output/snakemake_20260415_112218.log


Run this command to submit the job to the supercomputer:
  sbatch /users/9/user0011/ProjectName_AMF_output/submit_job.sh

```
Finally, if everything looks good, just run the `sbatch` command generated by your `run_dada2processing` command. This will submit your SLURM job to the cluster!

**A note about SBATCH/SLURM scripts if you are unfamiliar:** when you run `run_asv2otu`, it submits the job as a 'SLURM' submission to the computing cluster. This means the 'job' will get in a queue to eventually run. Depending on what time of day/week you submit it, it could take anywhere from 2 seconds to 30 minutes to initiate (usually towards the lower end in my experience). Once it starts, you will get an email saying it started. Next, you will either get an email that the script COMPLETED or FAILED. If it FAILED, it will specify an Exit Code number, usually 2. If it failed, contact me (Peter F) and I can help troubleshoot. Every time you run the pipeline, the job will output a `slurm_JOBID.out` and `slurm_JOBID.err` file into your project directory. You don't really need to worry about these files, UNLESS something fails, then both will be helpful for understanding what error was thrown/what went wrong.

After your job has started, you can watch its progress in real time using the following command:
```bash
# navigate to your project directory, then run:
tail -f slurm_*.out
```

This will show you a running output of what the pipeline is doing. It's not necessary to do this step, but if you want to track the pipeline in real time, this is an easy way to do it.

After the pipeline is done, or if you want to exit the live view of whats happening, hit Control+C & Enter on Mac (not sure what the equivalent is for Windows).


## Output Files:

The DADA2 pipeline will output a LOT of files, but you will likely be most interested in:

**For most datasets (ITS1, ITS2, 16S-V4, 18S-V4):**
1. *ProjectName__combined_sequences_ASVtaxa_DatabaseName.txt* - Your denoised ASV table with taxonomy assignments and bootstrap values
2. *ProjectName__OTUs_with_taxonomy_DatabaseName.txt* - Your final denoised OTU table with taxonomy assignments and boostrap values
2. *qc_summary.txt* - detailed QC file reporting number and perecentage of sequences removed in each step of the process

DADA2 also reports other QC statistics during its denoising process. To see the figures showing those denoising steps, check the *04_dada2_QCsummary/figures* folder. This contains 4 types of figures:

1. Error Model reporting for forward and reverse reads - the main thing you are looking for hear is that the black line mostly tracks the red line. Look online or in the *04_dada2_QCsummary/figures/README.md* for mor information
2. Quality Profile reporting for the first 3 read files (forward and reverse) - this shows the 'quality' of the read across it's length. Expect declines in quality at the end of the reverse read (don't ask me why, but it is a well-known phenomenon)
3. Read tracking: figure showing number of reads remaining at each denoising step.
4. Sequence length distribution: distribution of the length of your reads (hypothetically this should be unimodal and somewhat tightly clustered, however expect variation for ITS)


Each output folder (03_dada2, 04_dada2_QCsummary, 05_asv2otu) has associated READMEs within the folders which describe their contents. Please read those if you are interested in what is going on in the other output files.

## Summary of main pipeline steps:

### Step 1 - Adapter and Primer trimmming
To make sure that the reads you are processing only contain biological information, we need to make sure that adapters (used to communicate with the sequencing platform) and primers are removed from all sequences. In Illumina datasets, this is achieved by first searching and removing remaining adapters detected near the end of you forward and reverse reads using `cutadapt`, then primers are removed in a very similar fashion. In Aviti, the adapters are different (and also maybe change sometimes...), and adapter trimming has to be a little more stringent. First is an adapter trimming step using `trimmomatic`, then primers are trimmed using `cutadapt`, then adapters are search for again using `trimmomatic`. Trevor Gould says that he has found up to 4 or 5 copies of the Aviti adapters within reads, which is why this step is performed twice. He cautions however that sometimes reads can get chopped in half because primers are found deep within a read. Across both datasets, reads without primers detected are removed.

### Step 2 - DADA2 Denoising
DADA2 is platform designed for stepwise denoising of amplicon sequencing datasets (built to be optimized for Illumina, but we have adapted it for function with Aviti). The full details can be found in the original Callahan et al. paper, but the basic steps are as follows:
1. Filtering and Trimming: reads are filtered based on expected quality and length parameters specific to each primer set and sequencing platform
2. Dereplication: Reads are dereplicated based on exact identity and collapsed into single reads with counts.
3. Error Learning and Denoising: DADA2 reads your sequence files and learns the error patterns specific to that sequencing run. It then performs an additional denoising step based on this inference.
4. Merging: forward and reverse reads are merged, with minimum overlap requirements differing based on sequence dataset. NOTE: this is a step where a lot of reads can get thrown out because they cannot merge. This may happen because there is not enough overlap between the forward and reverse reads and is known to be a problem in 18S-V4 and 18S-AMF datasets. If you find that a low percentage of reads are merging, you may want to only use forward reads for your analysis.
5. Chimera Removal: Based on the error patterns learned by the DADA2 algorithm, it detects what it thinks are possible chimeric sequences based on error patterns and parent abundance. For complicated reasons related to Aviti error patterns, the minimum abundance of the parent required for a chimera merge/removal (minFoldParentOverAbundance) has been changed from 2 to 8 for Aviti datasets. Without these changes, a large proportion of real biological ASVs would be detected as chimeric incorrectly. Talk to Trevor Gould if you want more information on this.

### Step 3 - DADA2 RDP Taxonomy Assignment for ASVs
After the ASVs have been created and denoised, taxonomy is assigned using the *assignTaxonomy* function in DADA2, which uses the RDP Naive Bayesian classifier (doi: 10.1128/AEM.00062-07) to classify sequences. Taxonomy databases are unique to each primer (ITS1 and ITS2 have one database: UNITE sh) and are the same databases used by Trevor Gould to assign taxonomy. I have updated each taxonomy database to their most recent version, and they are located in the Kennedy lab shared taxonomy folder on the MSI computing cluster.

### Step 4 - ITSx (only relevant for ITS datasets)
 In essence, ITSx identifies and removes the highly conserved 5.8S part of the ITS2 barcodes. This makes for better clustering later, as the 5.8S regions inflates the similarity between sequences. It also discards sequences that are identified as not ITS, which may accidentally discard synthetic community members (synmock). If you run this and find that your synmock does not contain the expected high abundance synthetic community members, it may be because ITSx has removed them. Try re-running the pipeline without ITSx.

### Step 5 - VSEARCH Clustering

This step uses the VSEARCH algorithm to cluster ASVs at 97% similarity (weighted by sequence abundance; *see VSEARCH documentation*), then does an additional chimera check to find and remove any chimeras that might have slipped through DADA2. The number removed here should be quite small (check log file), but may vary wildly by primer set. 

**NOTE:** *Initial testing seems to suggest that up to 20% of the OTUs in 18S-V4 primer datasets are detected to be chimeric. My assumption with the 18S primer is that it is trying to pick up so much sequence diversity across the eukaryotic tree that it picks up a lot of errors a long the way. Nevertheless, I (PF) would like to do more testing to fully understand whats going on here.*

### Step 6 - MUMU/LULU OTU Curation

*mumu* is a post-clustering clean up algorithm that is supposed to find rare variants sequences of common sequences that leak through the clustering steps.
It works in two steps. All sequences are blasted against each other. Sequences with high similarity are checked for patterns in the OTU table. If a rare sequence with high similarity to a common sequence also shows a very similar occurrence pattern, it is merged with this "parent" sequence. These sequences are assumed to be additional sequencing errors that have made it through all previous cleaning steps. Mostly used to clean up singleton and doubleton sequences. *mumu* is a unix version of the LULU package (Frøslev et al., 2017), with some optimisation to run faster.

**NOTE:** *Similar to above, I have found in my preliminary testing of 18S-V4 datasets that more than 50% of OTUs are detected to have a parent OTU by mumu. Similar assumptions to above for why this is happening, but it is worth looking into this, and whether mumu parameters should be tweaked for 18S datasets.*

**Note 2:** *In preliminary testing, mumu was causing problematic merges of the ZyMock community in 16S-V4 bacterial datasets. This was due to a low blast % threshold considered for parent-daughter sequences, and a standard parent-daughter ratio of 1. To fix this, bacteria by default now run mumu with slightly different settings. The minimum BLAST % is 94% and the minimum ratio is 100. This is somewhat informed by the relatively more conserved nature of the 16S region, but is also somewhat arbitrarily chosen. Just expect that the mumu filtering step might merge comparatively less OTUs in 16S datasets compared to others. I do not expect this to change ecological inference. If you are still running into issues where the ZyMock community seems strange, contact me.*

### Step 7 - DADA2 RDP Taxonomy Assignment

After the OTUs have been created and curated, taxonomy is assigned to the centroid of each OTU using the exact same approach as above for ASV taxonomy assignment. 

This will likely take the longest of any step in the process -- somewhere between 5-45 minutes depending on how many OTUs there are, and the size of the reference database (PR2 takes the longest). 

**NOTE:** An important note about AMF taxonomy assignment - because the MaarjAM database includes ONLY AMF sequences, DADA2 can ONLY assign sequences as AMF. As such, a lot of OTUs which are definitely not AMF will be assigned to Mucoromycota/Glomeromycota with high bootstrap confidence. I would take all these assignments with a grain of salt, and my recommendation would be to only use Genus/Virtual Taxa assignments that have high confidence, or maybe Family-level assignments with high bootstrap. If you have suggestions for improving AMF filtering, please let me know as I am very open to suggestions!

## Citations

If you use this pipeline, please give credit to Luis Morgado for developing the original Zazzy Metabarcoding Pipeline (https://github.com/ekronold/Zazzy_metabarcoding_pipeline), and Eivind Kverme Ronold for making the initial updates for adaptation to DADA2 output ASV files. Please also give credit to Trevor Gould for testing and development of the adapter/primer trimming and DADA2 denoising settings, it would probably be helpful to cite his Aviti paper (currently in preprint) listed below.

*Additionally, please cite the primary packages and databases used in the pipeline:*

### Methods:
Gould, T. J., Taylor, M., & Santelli, C. (2026). Aviti Sequencing and Marker Gene Data Analysis (p. 2026.02.06.704475). bioRxiv. https://doi.org/10.64898/2026.02.06.704475


### Packages:

ITSx: Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for use in environmental sequencing. Johan Bengtsson-Palme, Vilmar Veldre, Martin Ryberg, Martin Hartmann, Sara Branco, Zheng Wang, Anna Godhe, Yann Bertrand, Pierre De Wit, Marisol Sanchez, Ingo Ebersberger, Kemal Sanli, Filipe de Souza, Erik Kristiansson, Kessy Abarenkov, K. Martin Eriksson, R. Henrik Nilsson *Methods in Ecology and Evolution, 4: 914-919, 2013* (https://doi.org/10.1111/2041-210X.12073)

Rognes T, Flouri T, Nichols B, Quince C, Mahé F. (2016) VSEARCH: a versatile open source tool for metagenomics. PeerJ 4:e2584.
doi: [10.7717/peerj.2584](https://doi.org/10.7717/peerj.2584)

Frédéric Mahé. (2026) mumu: post-clustering curation tool for metabarcoding data v.1.1.2 (https://github.com/frederic-mahe/mumu)

Frøslev, T. G., Kjøller, R., Bruun, H. H., Ejrnæs, R., Brunbjerg, A. K., Pietroni, C., & Hansen, A. J. (2017). Algorithm for post-clustering curation of DNA amplicon data yields reliable biodiversity estimates. Nature Communications, 8(1), 1188. (https://doi.org/10.1038/s41467-017-01312-x)

Callahan, B., McMurdie, P., Rosen, M. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016). https://doi.org/10.1038/nmeth.3869

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1), 10–12. doi:10.14806/ej.17.1.200

Anthony M. Bolger, Marc Lohse, Bjoern Usadel, Trimmomatic: a flexible trimmer for Illumina sequence data, Bioinformatics, Volume 30, Issue 15, August 2014, Pages 2114–2120, https://doi-org.ezp1.lib.umn.edu/10.1093/bioinformatics/btu170


### Databases:

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C. et al. 2013. The Protist Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Res. 41:D597–604.

Chuvochina M, Gerken J, Frentrup M, Sandikci Y, Goldmann R, Freese HM, Göker M, Sikorski J, Yarza P, Quast C, Peplies J, Glöckner FO, Reimer LC (2026) SILVA in 2026: a global core biodata resource for rRNA within the DSMZ digital diversity. Nucleic Acids Research, gkaf1247.

Abarenkov K, Nilsson RH, Larsson K-H, Taylor AFS, May TW, Frøslev TG, Pawlowska J, Lindahl B, Põldmaa K, Truong C, Vu D, Hosoya T, Niskanen T, Piirmann T, Ivanov F, Zirk A, Peterson M, Cheeke TE, Ishigami Y, Jansson AT, Jeppesen TS, Kristiansson E, Mikryukov V, Miller JT, Oono R, Ossandon FJ, Paupério J, Saar I, Schigel D, Suija A, Tedersoo L, Kõljalg U. 2024. The UNITE database for molecular identification and taxonomic communication of fungi and other eukaryotes: sequences, taxa and classifications reconsidered. Nucleic Acids Research, https://doi.org/10.1093/nar/gkad1039

Öpik, M., Vanatoa, A., Vanatoa, E., Moora, M., Davison, J., Kalwij, J.M., Reier, Ü., Zobel, M. 2010. The online database MaarjAM reveals global and ecosystemic distribution patterns in arbuscular mycorrhizal fungi (Glomeromycota). New Phytologist 188: 223-241.

Leho Tedersoo, Mahdieh S Hosseyni Moghaddam, Vladimir Mikryukov, Ali Hakimzadeh, Mohammad Bahram, R Henrik Nilsson, Iryna Yatsiuk, Stefan Geisen, Arne Schwelm, Kasia Piwosz, Marko Prous, Sirje Sildever, Dominika Chmolowska, Sonja Rueckert, Pavel Skaloud, Peeter Laas, Marco Tines, Jae-Ho Jung, Ji Hye Choi, Saad Alkahtani, Sten Anslan , EUKARYOME: the rRNA gene reference database for identification of all eukaryotes, Database, Volume 2024, 2024, baae043, https://doi.org/10.1093/database/baae043

## License

Internal use only
