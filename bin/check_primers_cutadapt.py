#!/usr/bin/env python3
"""
Primer detection using cutadapt
Searches raw FASTQ files for all known primers and reports which ones are found
Usage: python3 check_primers_cutadapt.py <fastq_dir> <expected_fwd_primer> <expected_rev_primer> <amplicon>
"""

import subprocess
import json
import sys
import os
import tempfile
import glob
from pathlib import Path

# Primer reference database
PRIMER_FWD = {
    "ITS1F": "CTTGGTCATTTAGAGGAAGTAA",
    "ITS1F_KYO1": "CTHGGTCATTTAGAGGAASTAA",
    "5.8S-FUN": "AACTTTYRRCAAYGGATCWCT",
    "fITS7_1": "GTGARTCATCGAATCTTG",
    "5.8SR": "TCGATGAAGAACGCAGCG",
    "WANDA_1": "CAGCCGCGGTAATTCCAGCT",
    "515f": "GTGYCAGCMGCCGCGGTAA",
    "27F": "AGRTTTGATYMTGGCTCAG",
    "F04": "GCTTGTCTCAAAGATTAAGCC",
    "616f": "TTAAARVGYTCGTAGTYG",
}

PRIMER_REV = {
    "ITS2": "GCTGCGTTCTTCATCGATGC",
    "ITS2_KYO2": "TTYRCTRYGTTCTTCATC",
    "ITS4-Fun": "AGCCTCCGCTTATTGATATGCTTAART",
    "ITS4_1": "TCCTCCGCTTATTGATATGC",
    "ITS4": "TCCTCCGCTTATTGATATGC",
    "AML2_1": "GAACCCAAACACTTTGGTTTCC",
    "806r": "GGACTACNVGGGTWTCTAAT",
    "1492R": "RGYTACCTTGTTACGACTT",
    "R22": "GCCTGCTGCCTTCCTTGGA",
    "1132r": "CCGTCAATTHCTTTYAART",
}


def find_fastq_files(fastq_dir):
    """Find first non-NC R1 and R2 files"""
    r1_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
    r2_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

    # Filter out negative controls (contain NC in the filename before first underscore)
    # Matches patterns like: FAB2NC-V4_S95, SampleNC_S1, etc.
    r1_files = [f for f in r1_files if "NC" not in os.path.basename(f).split("_")[0]]
    r2_files = [f for f in r2_files if "NC" not in os.path.basename(f).split("_")[0]]

    if not r1_files or not r2_files:
        return None, None

    return r1_files[0], r2_files[0]


def run_cutadapt_primer_search(fastq_file, primers_dict, is_reverse=False):
    """
    Run cutadapt to search for primers in FASTQ file
    Returns dict of {primer_name: match_count}
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write primers to FASTA
        fasta_file = os.path.join(tmpdir, "primers.fasta")
        with open(fasta_file, "w") as f:
            for name, seq in primers_dict.items():
                f.write(f">{name}\n{seq}\n")

        json_output = os.path.join(tmpdir, "results.json")

        # Run cutadapt
        # -g for forward strand (5' adapters)
        # -G for reverse strand (3' adapters on reverse read)
        flag = "-G" if is_reverse else "-g"

        cmd = [
            "cutadapt",
            "--no-trim",
            "--untrimmed-output=/dev/null",
            f"{flag}",
            f"file:{fasta_file}",
            f"--json={json_output}",
            fastq_file,
        ]

        try:
            subprocess.run(cmd, capture_output=True, check=True, timeout=300)
        except subprocess.CalledProcessError as e:
            print(f"ERROR: cutadapt failed: {e.stderr.decode()}", file=sys.stderr)
            sys.exit(1)

        # Parse JSON output
        with open(json_output) as f:
            data = json.load(f)

        # Debug: print JSON structure if read_counts is missing
        if "read_counts" not in data:
            print(f"DEBUG: JSON structure: {json.dumps(data, indent=2)[:500]}", file=sys.stderr)

        # Extract adapter matches
        results = {}

        # Handle different cutadapt JSON structures
        # read_counts can be either:
        # - list: [input, output]
        # - dict: {"input": N, "output": N}
        # - direct int in older versions
        read_counts_data = data.get("read_counts")
        if isinstance(read_counts_data, list) and len(read_counts_data) > 0:
            total_reads = read_counts_data[0]
        elif isinstance(read_counts_data, dict):
            total_reads = read_counts_data.get("input", read_counts_data.get("total", 0))
        elif isinstance(read_counts_data, int):
            total_reads = read_counts_data
        else:
            # Fallback: count from adapters if present
            total_reads = sum(a.get("matches", 0) for a in data.get("adapters", [])) or 1

        for adapter_info in data.get("adapters", []):
            name = adapter_info["name"]
            matches = adapter_info["matches"]
            results[name] = {"matches": matches, "pct": 100 * matches / total_reads if total_reads > 0 else 0}

        return results, total_reads


def main():
    if len(sys.argv) < 5:
        print("Usage: check_primers_cutadapt.py <fastq_dir> <expected_fwd> <expected_rev> <amplicon>")
        sys.exit(1)

    fastq_dir = sys.argv[1]
    expected_fwd = sys.argv[2]
    expected_rev = sys.argv[3]
    amplicon = sys.argv[4]

    print("Primer Detection with Cutadapt")
    print("=" * 50)
    print(f"Fastq directory: {fastq_dir}")
    print(f"Amplicon:        {amplicon}")
    print(f"Expected fwd:    {expected_fwd}")
    print(f"Expected rev:    {expected_rev}")
    print()

    # Find files
    r1_file, r2_file = find_fastq_files(fastq_dir)
    if not r1_file or not r2_file:
        print("ERROR: Could not find non-NC R1/R2 files in", fastq_dir)
        sys.exit(1)

    print("Checking files (excluding negative controls):")
    print(f"  R1: {os.path.basename(r1_file)}")
    print(f"  R2: {os.path.basename(r2_file)}")
    print()

    # Forward primer analysis
    print("FORWARD PRIMER ANALYSIS")
    print("-" * 50)
    print(f"Expected: {expected_fwd}")

    fwd_results, fwd_total = run_cutadapt_primer_search(r1_file, PRIMER_FWD, is_reverse=False)

    # Sort by match frequency
    sorted_fwd = sorted(fwd_results.items(), key=lambda x: x[1]["matches"], reverse=True)

    print(f"Total reads processed: {fwd_total}")
    print("Results (sorted by frequency):")

    for name, data in sorted_fwd:
        matches = data["matches"]
        pct = data["pct"]
        marker = ""

        # Mark expected primer
        if name == "515f" and expected_fwd == PRIMER_FWD.get("515f"):
            marker = " <-- EXPECTED"
        elif name in PRIMER_FWD and PRIMER_FWD[name] == expected_fwd:
            marker = " <-- EXPECTED"

        # Highlight unexpected high-frequency primers
        if pct > 80 and not marker:
            marker = " <-- UNEXPECTED! (>80%)"
        elif pct > 50 and not marker:
            marker = " <-- unexpected (>50%)"

        print(f"  {name:15} {matches:7} reads ({pct:5.1f}%){marker}")

    # Reverse primer analysis
    print()
    print("REVERSE PRIMER ANALYSIS")
    print("-" * 50)
    print(f"Expected: {expected_rev}")

    rev_results, rev_total = run_cutadapt_primer_search(r2_file, PRIMER_REV, is_reverse=True)

    # Sort by match frequency
    sorted_rev = sorted(rev_results.items(), key=lambda x: x[1]["matches"], reverse=True)

    print(f"Total reads processed: {rev_total}")
    print("Results (sorted by frequency):")

    for name, data in sorted_rev:
        matches = data["matches"]
        pct = data["pct"]
        marker = ""

        # Mark expected primer
        if name == "806r" and expected_rev == PRIMER_REV.get("806r"):
            marker = " <-- EXPECTED"
        elif name in PRIMER_REV and PRIMER_REV[name] == expected_rev:
            marker = " <-- EXPECTED"

        # Highlight unexpected high-frequency primers
        if pct > 80 and not marker:
            marker = " <-- UNEXPECTED! (>80%)"
        elif pct > 50 and not marker:
            marker = " <-- unexpected (>50%)"

        print(f"  {name:15} {matches:7} reads ({pct:5.1f}%){marker}")

    # Summary
    print()
    print("SUMMARY & RECOMMENDATIONS")
    print("=" * 50)

    # Find expected primers in results
    expected_fwd_found = any(
        PRIMER_FWD.get(name) == expected_fwd and data["pct"] > 80 for name, data in fwd_results.items()
    )
    expected_rev_found = any(
        PRIMER_REV.get(name) == expected_rev and data["pct"] > 80 for name, data in rev_results.items()
    )

    if expected_fwd_found and expected_rev_found:
        print("✓ Both expected primers detected at >80% frequency. Good!")
    elif expected_fwd_found or expected_rev_found:
        print("⚠ One expected primer found at high frequency.")
        print("  Review the other direction for issues.")
    else:
        print("✗ Expected primers NOT detected at high frequency.")
        print("  REVIEW RESULTS ABOVE for alternative primers detected.")

    # Check for unexpected high-frequency primers
    unexpected_fwd = [
        (name, data)
        for name, data in fwd_results.items()
        if data["pct"] > 80 and (not any(PRIMER_FWD.get(n) == expected_fwd for n in [name]))
    ]
    unexpected_rev = [
        (name, data)
        for name, data in rev_results.items()
        if data["pct"] > 80 and (not any(PRIMER_REV.get(n) == expected_rev for n in [name]))
    ]

    if unexpected_fwd:
        print()
        print("⚠ Unexpected forward primers detected at >80%:")
        for name, data in unexpected_fwd:
            print(f"  - {name}: {data['pct']:.1f}%")

    if unexpected_rev:
        print()
        print("⚠ Unexpected reverse primers detected at >80%:")
        for name, data in unexpected_rev:
            print(f"  - {name}: {data['pct']:.1f}%")

    print()
    print("Do NOT proceed with SLURM submission until primers are confirmed!")


if __name__ == "__main__":
    main()
