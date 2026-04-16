#!/usr/bin/env python3
"""
Adapter auto-detection using cutadapt's --detect-adapters
Identifies adapter sequences present in FASTQ files
Usage: python3 detect_adapters.py <fastq_dir> <num_samples> <output_dir>
"""

import subprocess
import json
import sys
import os
import glob
import tempfile
from pathlib import Path

# Default Nextera XT adapters (fallback if detection fails)
DEFAULT_NEXTERA_ADAPTERS = [
    "ATCTCGTATGCCGTCTTCTGCTTG",      # i7 Nextera XT
    "CAAGCAGAAGACGGCATACGAGAT",      # i7 reverse complement
    "GTGTAGATCTCGGTGGTCGCCGTATCATT", # i5 Nextera XT
    "AATGATACGGCGACCACCGAGATCTACAC", # i5 reverse complement
]


def find_fastq_files(fastq_dir, num_samples=3):
    """Find first N non-NC R1 and R2 file pairs"""
    r1_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R1_001.fastq.gz")))
    r2_files = sorted(glob.glob(os.path.join(fastq_dir, "*_R2_001.fastq.gz")))

    # Filter out negative controls (contain NC in the filename before first underscore)
    r1_files = [f for f in r1_files if "NC" not in os.path.basename(f).split("_")[0]]
    r2_files = [f for f in r2_files if "NC" not in os.path.basename(f).split("_")[0]]

    if not r1_files or not r2_files:
        return []

    # Match R1 and R2 files by sample name stem (not by sorted index)
    r1_stems = {os.path.basename(f).replace("_R1_001.fastq.gz", ""): f for f in r1_files}
    r2_stems = {os.path.basename(f).replace("_R2_001.fastq.gz", ""): f for f in r2_files}
    matched_stems = sorted(set(r1_stems) & set(r2_stems))

    # Return first N matched pairs
    pairs = []
    for stem in matched_stems[:num_samples]:
        pairs.append((r1_stems[stem], r2_stems[stem]))

    return pairs


def run_detect_adapters(r1_file, r2_file, read_limit=100000):
    """
    Run cutadapt --detect-adapters on paired-end reads
    Returns dict with detected adapters and statistics
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        json_output = os.path.join(tmpdir, "results.json")

        # Run cutadapt with --detect-adapters
        # Process first N reads to speed up detection
        cmd = [
            "cutadapt",
            "--detect-adapters",
            "--max-reads", str(read_limit),
            "--json", json_output,
            "-o", "/dev/null",
            "-p", "/dev/null",
            r1_file,
            r2_file,
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, check=True, timeout=300, text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"ERROR: cutadapt --detect-adapters failed: {e.stderr}", file=sys.stderr)
            return None
        except subprocess.TimeoutExpired:
            print(f"ERROR: cutadapt --detect-adapters timed out", file=sys.stderr)
            return None

        # Parse JSON output
        try:
            with open(json_output) as f:
                data = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            print(f"ERROR: Failed to parse cutadapt JSON: {e}", file=sys.stderr)
            return None

        # Extract detected adapters (cutadapt JSON uses adapters_read1/adapters_read2)
        result_data = {
            "r1_adapters": data.get("adapters_read1", []),
            "r2_adapters": data.get("adapters_read2", []),
        }

        return result_data


def main():
    if len(sys.argv) < 4:
        print("Usage: detect_adapters.py <fastq_dir> <num_samples> <output_dir>")
        sys.exit(1)

    fastq_dir = sys.argv[1]
    num_samples = int(sys.argv[2])
    output_dir = sys.argv[3]

    print("Adapter Auto-Detection with Cutadapt")
    print("=" * 60)
    print(f"Fastq directory: {fastq_dir}")
    print(f"Samples to scan: {num_samples}")
    print()

    # Find files
    sample_pairs = find_fastq_files(fastq_dir, num_samples)
    if not sample_pairs:
        print("ERROR: Could not find non-NC R1/R2 files in", fastq_dir)
        sys.exit(1)

    print(f"Found {len(sample_pairs)} non-NC sample pair(s) to scan:")
    for r1, r2 in sample_pairs:
        print(f"  - {os.path.basename(r1)}")
    print()

    # Run detection on each sample pair
    all_r1_adapters = []
    all_r2_adapters = []

    for i, (r1_file, r2_file) in enumerate(sample_pairs, 1):
        sample_name = os.path.basename(r1_file).replace("_R1_001.fastq.gz", "")
        print(f"[Sample {i}/{len(sample_pairs)}] Detecting adapters in {sample_name}...")

        result = run_detect_adapters(r1_file, r2_file)
        if result is None:
            print(f"  ⚠ Detection failed for {sample_name}")
            continue

        r1_adapters = result.get("r1_adapters", [])
        r2_adapters = result.get("r2_adapters", [])

        if r1_adapters:
            print(f"  R1 adapters detected: {len(r1_adapters)}")
            for adapter in r1_adapters:
                print(f"    - {adapter}")
            all_r1_adapters.extend(r1_adapters)

        if r2_adapters:
            print(f"  R2 adapters detected: {len(r2_adapters)}")
            for adapter in r2_adapters:
                print(f"    - {adapter}")
            all_r2_adapters.extend(r2_adapters)

        if not r1_adapters and not r2_adapters:
            print(f"  ℹ No adapters detected in {sample_name}")

        print()

    # Deduplicate detected adapters while preserving order
    unique_r1 = []
    seen = set()
    for adapter in all_r1_adapters:
        if adapter not in seen:
            unique_r1.append(adapter)
            seen.add(adapter)

    unique_r2 = []
    seen = set()
    for adapter in all_r2_adapters:
        if adapter not in seen:
            unique_r2.append(adapter)
            seen.add(adapter)

    # Summary and output
    print("=" * 60)
    print("ADAPTER DETECTION SUMMARY")
    print("=" * 60)

    if unique_r1 or unique_r2:
        print(f"✓ Detected adapters across {len(sample_pairs)} sample(s):")
        print()

        if unique_r1:
            print(f"R1 adapters ({len(unique_r1)}):")
            for adapter in unique_r1:
                print(f"  -a {adapter}")
            print()

        if unique_r2:
            print(f"R2 adapters ({len(unique_r2)}):")
            for adapter in unique_r2:
                print(f"  -A {adapter}")
            print()

        print("These WILL BE COMBINED with default Nextera XT adapters:")
        print("Default adapters being used:")
        for adapter in DEFAULT_NEXTERA_ADAPTERS:
            print(f"  - {adapter}")
        print()

    else:
        print("ℹ No adapters detected via --detect-adapters")
        print("Will use default Nextera XT adapters:")
        for adapter in DEFAULT_NEXTERA_ADAPTERS:
            print(f"  - {adapter}")
        print()

    # Build final adapter list (detected + defaults)
    final_adapters = DEFAULT_NEXTERA_ADAPTERS.copy()
    for adapter in unique_r1 + unique_r2:
        if adapter not in final_adapters:
            final_adapters.append(adapter)

    # Save to JSON for run_dada2processing to read
    adapter_detection = {
        "detected_r1": unique_r1,
        "detected_r2": unique_r2,
        "default_adapters": DEFAULT_NEXTERA_ADAPTERS,
        "final_adapters": final_adapters,
        "num_samples_scanned": len(sample_pairs),
    }

    json_output_path = os.path.join(output_dir, ".adapter_detection.json")
    with open(json_output_path, "w") as f:
        json.dump(adapter_detection, f, indent=2)

    print(f"[Adapter detection results saved to {json_output_path}]")


if __name__ == "__main__":
    main()
