#!/usr/bin/env python3
"""
Extract primer sequences from primers.json config.
Used by bin/run_dada2processing to replace broken inline Python heredoc.

Usage:
    python3 get_primers.py --primers-json <path> --amplicon <name>

Output:
    Prints two space-separated sequences: <fwd_seq> <rev_seq>
    Exits nonzero if amplicon not found.
"""

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser(description="Extract primer sequences from JSON config")
    parser.add_argument("--primers-json", required=True, help="Path to primers.json")
    parser.add_argument("--amplicon", required=True, help="Amplicon name (e.g., 16S-V4)")
    args = parser.parse_args()

    try:
        with open(args.primers_json) as f:
            primers = json.load(f)
    except FileNotFoundError:
        print(f"ERROR: primers file not found: {args.primers_json}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"ERROR: invalid JSON in {args.primers_json}: {e}", file=sys.stderr)
        sys.exit(1)

    if args.amplicon not in primers:
        available = ", ".join(sorted(primers.keys()))
        print(f"ERROR: amplicon '{args.amplicon}' not found in {args.primers_json}", file=sys.stderr)
        print(f"Available amplicons: {available}", file=sys.stderr)
        sys.exit(1)

    amplicon_config = primers[args.amplicon]
    if "fwd" not in amplicon_config or "rev" not in amplicon_config:
        print(f"ERROR: amplicon '{args.amplicon}' missing 'fwd' or 'rev' key", file=sys.stderr)
        sys.exit(1)

    fwd = amplicon_config["fwd"]
    rev = amplicon_config["rev"]
    print(f"{fwd} {rev}")


if __name__ == "__main__":
    main()
