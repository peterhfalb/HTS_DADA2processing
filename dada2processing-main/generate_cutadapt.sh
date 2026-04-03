#!/bin/bash

# Get region name from the first command-line argument
REGION="$1"

case "$REGION" in
  ITS1)
    # ITS1
    for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed -a ^CTTGGTCATTTAGAGGAAGTAA...GCATCGATGAAGAACGCAGC -A ^GCTGCGTTCTTCATCGATGC...TTACTTCCTCTAAATGACCAAG -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  ITS2)
    # ITS2
    for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed  -a ^TCGATGAAGAACGCAGCG...GCATATCAATAAGCGGAGGA -A ^TCCTCCGCTTATTGATATGC...CGCTGCGTTCTTCATCGA -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  V4)
    # V4
	for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed  -a ^GTGCCAGCNGCCGCGGTAA...ATTAGANACCCNNGTAGTCC -A ^GGACTACNNGGGTNTCTAAT...TTACCGCGGCNGCTGGCAC -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  V3V5)
    # V3V5
	for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed  -a ^CCTACGGGNGGCNGCAG...AAACTNAAANNAATTGNCGG -A ^CCGNCAATTNNTTTNAGTTT...CTGCNGCCNCCCGTAGG -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  V5V6)
    # V5V6
	for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed  -a ^NGGATTAGATACCC...AGGTGNTGCATGGNNGTCG -A ^CGACNNCCATGCANCACCT...GGGTATCTAATCCN -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  V3V4)
    # V3V4
	for i in *_R1_001.paired.fastq.gz; do echo "cutadapt --cores 8 --pair-filter=any --minimum-length 100 --discard-untrimmed  -a ^CCTACGGGNGGCNGCAG...ATTAGANACCCNNGTAGTCC -A ^GGACTACNNGGGTNTCTAAT...CTGCNGCCNCCCGTAGG -o ../02_filtered/${i//_R1_001.paired.fastq.gz/_R1_001.fastq.gz} -p ../02_filtered/${i//_R1_001.paired.fastq.gz/_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.paired.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.cmd; done
    ;;
  *)
    echo "Error: REGION must be one of: ITS1, ITS2, V4, V3V5, V5V6, V3V4" >&2
    exit 1
    ;;
esac
