#!/usr/bin/env python3
"""
Aggregate *.summary.tsv files from QC workflow runs into a single TSV.

Finds all *.summary.tsv files under the current directory (or a given path),
writes the header once, then appends one data row per file.

Usage:
    python aggregate_qc_summary.py
    python aggregate_qc_summary.py --search_dir /path/to/batch
    python aggregate_qc_summary.py --output my_batch_summary.tsv
"""

import argparse
import glob
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--search_dir", default=".",
        help="Directory to search for *.summary.tsv files (default: current directory)"
    )
    parser.add_argument(
        "--output", default="summary_agg.tsv",
        help="Output file name (default: summary_agg.tsv)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    files = sorted(glob.glob(os.path.join(args.search_dir, "**", "*.summary.tsv"), recursive=True))

    if not files:
        sys.exit(f"No *.summary.tsv files found under {args.search_dir}")

    with open(args.output, "w") as out:
        for i, path in enumerate(files):
            with open(path) as fh:
                lines = fh.read().splitlines()
            if not lines:
                print(f"Warning: empty file skipped: {path}")
                continue
            if i == 0:
                out.write(lines[0] + "\n")  # header from first file
            if len(lines) > 1:
                out.write(lines[1] + "\n")  # data row
            else:
                print(f"Warning: no data row in {path}")

    print(f"Aggregated {len(files)} file(s) into {args.output}")


if __name__ == "__main__":
    main()
