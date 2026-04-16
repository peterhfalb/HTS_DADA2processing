#!/usr/bin/env python3
"""
Primer detection using cutadapt
Searches raw FASTQ files for all known primers and reports which ones are found
Usage: python3 check_primers_cutadapt.py <fastq_dir> <expected_fwd_primer> <expected_rev_primer> <amplicon> <output_dir>
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


def run_cutadapt_primer_search(r1_file, r2_file, primers_dict, is_reverse=False):
    """
    Run cutadapt to search for primers in paired-end FASTQ files
    Returns dict of {primer_name: match_count}
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write primers to FASTA
        fasta_file = os.path.join(tmpdir, "primers.fasta")
        with open(fasta_file, "w") as f:
            for name, seq in primers_dict.items():
                f.write(f">{name}\n{seq}\n")

        json_output = os.path.join(tmpdir, "results.json")

        # Run cutadapt in paired-end mode
        # -g for forward strand (5' adapters on R1)
        # -G for reverse strand (3' adapters on R2)
        flag = "-G" if is_reverse else "-g"

        cmd = [
            "cutadapt",
            "--no-trim",
            f"{flag}",
            f"file:{fasta_file}",
            f"--json={json_output}",
            "-o", "/dev/null",  # Discard trimmed R1 output
            "-p", "/dev/null",  # Discard trimmed R2 output
            r1_file,
            r2_file,
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
            total_reads = 1

        # Parse adapters from either old format (adapters) or new format (adapters_read1/adapters_read2)
        adapters_list = []
        if "adapters" in data:
            adapters_list = data.get("adapters", [])
        elif is_reverse and "adapters_read2" in data:
            adapters_list = data.get("adapters_read2", [])
        elif not is_reverse and "adapters_read1" in data:
            adapters_list = data.get("adapters_read1", [])

        for adapter_info in adapters_list:
            name = adapter_info["name"]
            # Handle both old "matches" and new "total_matches" field names
            matches = adapter_info.get("total_matches") or adapter_info.get("matches", 0)
            results[name] = {"matches": matches, "pct": 100 * matches / total_reads if total_reads > 0 else 0}

        return results, total_reads


def main():
    if len(sys.argv) < 6:
        print("Usage: check_primers_cutadapt.py <fastq_dir> <expected_fwd> <expected_rev> <amplicon> <output_dir>")
        sys.exit(1)

    fastq_dir = sys.argv[1]
    expected_fwd = sys.argv[2]
    expected_rev = sys.argv[3]
    amplicon = sys.argv[4]
    output_dir = sys.argv[5]

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

    fwd_results, fwd_total = run_cutadapt_primer_search(r1_file, r2_file, PRIMER_FWD, is_reverse=False)

    # Sort by match frequency
    sorted_fwd = sorted(fwd_results.items(), key=lambda x: x[1]["matches"], reverse=True)

    print(f"Total reads processed: {fwd_total}")
    print("Results (sorted by frequency):")

    for name, data in sorted_fwd:
        matches = data["matches"]
        pct = data["pct"]
        marker = ""

        # Mark expected primer (exact match OR primer contains expected sequence)
        if name in PRIMER_FWD:
            if PRIMER_FWD[name] == expected_fwd:
                marker = " <-- EXPECTED"
            elif expected_fwd in PRIMER_FWD[name]:
                # Detected primer is a superset/degenerate version of expected
                marker = " <-- EXPECTED (degenerate variant)"

        # Highlight unexpected high-frequency primers (but not if they contain expected sequence)
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

    rev_results, rev_total = run_cutadapt_primer_search(r1_file, r2_file, PRIMER_REV, is_reverse=True)

    # Sort by match frequency
    sorted_rev = sorted(rev_results.items(), key=lambda x: x[1]["matches"], reverse=True)

    print(f"Total reads processed: {rev_total}")
    print("Results (sorted by frequency):")

    for name, data in sorted_rev:
        matches = data["matches"]
        pct = data["pct"]
        marker = ""

        # Mark expected primer (exact match OR primer contains expected sequence)
        if name in PRIMER_REV:
            if PRIMER_REV[name] == expected_rev:
                marker = " <-- EXPECTED"
            elif expected_rev in PRIMER_REV[name]:
                # Detected primer is a superset/degenerate version of expected
                marker = " <-- EXPECTED (degenerate variant)"

        # Highlight unexpected high-frequency primers (but not if they contain expected sequence)
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
    # Check for exact match OR if expected primer is contained in detected primer (for degenerate primers)
    expected_fwd_found = any(
        (PRIMER_FWD.get(name) == expected_fwd or expected_fwd in PRIMER_FWD.get(name, "")) and data["pct"] > 80
        for name, data in fwd_results.items()
    )
    expected_rev_found = any(
        (PRIMER_REV.get(name) == expected_rev or expected_rev in PRIMER_REV.get(name, "")) and data["pct"] > 80
        for name, data in rev_results.items()
    )

    # Special note for 18S-V4
    if amplicon == "18S-V4":
        print()
        print("⚠️  NOTE: 18S-V4 DATASETS")
        print("-" * 50)
        print("18S-V4 reverse reads (R2) are known to have quality issues:")
        print("  - Lower quality compared to forward reads")
        print("  - Inconsistent reverse primer binding")
        print("  - Variable primer trimming efficiency")
        print()
        print("If reverse primer matching is <80% or scattered across multiple")
        print("primers, this is expected behavior and does NOT indicate a problem.")
        print()
        print("⚠️  IMPORTANT: Check the merge rate in QC stats!")
        print("If merge rate is very poor (<50%), consider using forward reads only")
        print("(use --fwd-reads-only flag in run_dada2processing command) rather than attempting to merge.")
        print()

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

    # Output detected primers as JSON for automated primer switching
    # Find the highest-frequency primer for each direction
    detected_primers = {
        "fwd_name": sorted_fwd[0][0] if sorted_fwd else None,
        "fwd_seq": PRIMER_FWD.get(sorted_fwd[0][0]) if sorted_fwd else None,
        "fwd_pct": sorted_fwd[0][1]["pct"] if sorted_fwd else 0,
        "rev_name": sorted_rev[0][0] if sorted_rev else None,
        "rev_seq": PRIMER_REV.get(sorted_rev[0][0]) if sorted_rev else None,
        "rev_pct": sorted_rev[0][1]["pct"] if sorted_rev else 0,
        "expected_fwd": expected_fwd,
        "expected_rev": expected_rev,
    }

    # Write detected primers to a JSON file for the bash script to read
    # Save to output_dir to avoid permission issues with shared FASTQ directories
    json_output_path = os.path.join(output_dir, ".primer_detection.json")
    with open(json_output_path, "w") as f:
        json.dump(detected_primers, f, indent=2)

    print(f"\n[Primer detection results saved]")


if __name__ == "__main__":
    main()
