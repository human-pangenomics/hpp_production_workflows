#!/usr/bin/env python3
"""
Generate input JSONs and samples CSV for hic_qc_workflow, then print the sbatch command.

Usage:
    python setup_hic_qc.py files.txt
    python setup_hic_qc.py files.txt --job_name my_batch --partition long

Input file format (tab-separated R1 and R2, one lane per line):
    /path/to/sample1_R1.fastq.gz    /path/to/sample1_R2.fastq.gz
    /path/to/sample2_R1.fastq.gz    /path/to/sample2_R2.fastq.gz

The file_id is derived from the R1 filename.
"""

import argparse
import json
import os
import re
import csv
import sys


WDL = "/private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/hic_qc_workflow.wdl"
TOIL_SCRIPT = "/private/nanopore/hprc_qc/scripts/toil_sbatch_slurm.sh"

HPRC_ID_RE = re.compile(r'((?:HG|NA|PG)\d{5})')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "file_list",
        help="Tab-separated text file with one lane per line: r1_path<TAB>r2_path"
    )
    parser.add_argument(
        "--job_name", default="hic_qc",
        help="Slurm job name (default: hic_qc)"
    )
    parser.add_argument(
        "--partition", default="long",
        help="Slurm partition (default: long)"
    )
    parser.add_argument(
        "--wdl", default=WDL,
        help=f"Path to hic_qc_workflow.wdl (default: {WDL})"
    )
    return parser.parse_args()


def strip_ext(path):
    name = os.path.basename(path)
    for ext in (".fastq.gz", ".fastq", ".fq.gz", ".fq", ".bam"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def normalize_sample_id(hprc_id):
    """Normalize non-standard prefixes to HG or NA. PG -> HG (same individuals)."""
    if hprc_id.startswith("PG"):
        return "HG" + hprc_id[2:]
    return hprc_id


def build_sample_id(path):
    match = HPRC_ID_RE.search(path)
    basename = strip_ext(path)
    if match:
        hprc_id = normalize_sample_id(match.group(1))
        if basename.startswith(hprc_id):
            return basename
        return f"{hprc_id}_{basename}"
    print(f"Warning: could not parse HPRC sample ID from path, using filename: {os.path.basename(path)}", file=sys.stderr)
    return basename


def main():
    args = parse_args()

    with open(args.file_list) as fh:
        lines = [line.strip() for line in fh if line.strip()]

    if not lines:
        sys.exit(f"No entries found in {args.file_list}")

    samples = []
    for i, line in enumerate(lines, 1):
        parts = line.split("\t")
        if len(parts) != 2:
            sys.exit(f"Line {i} must have exactly two tab-separated fields (R1<TAB>R2):\n  {line}")
        r1, r2 = parts[0].strip(), parts[1].strip()
        samples.append((r1, r2))

    cwd = os.getcwd()
    json_dir = os.path.join(cwd, "input_jsons")
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(os.path.join(cwd, "slurm_logs"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "analysis"), exist_ok=True)

    csv_rows = []
    for r1, r2 in samples:
        sid = build_sample_id(r1)
        entry = {
            "hic_qc_wf.hic_reads": [r1, r2],
            "hic_qc_wf.file_id": sid,
        }
        json_path = os.path.join(json_dir, f"{sid}_hic_qc_workflow.json")
        with open(json_path, "w") as fh:
            json.dump(entry, fh, indent=2)
        csv_rows.append({"sample_id": sid, "r1": r1, "r2": r2})

    csv_path = os.path.join(cwd, "samples.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "r1", "r2"])
        writer.writeheader()
        writer.writerows(csv_rows)

    n = len(samples)
    print(f"\nCreated {n} input JSON(s) in {json_dir}/")
    print(f"Created {csv_path}\n")
    print("Run the following to submit:\n")
    print(
        f"sbatch \\\n"
        f"    --job-name={args.job_name} \\\n"
        f"    --array=[1-{n}]%{n} \\\n"
        f"    --partition={args.partition} \\\n"
        f"    {TOIL_SCRIPT} \\\n"
        f"    --wdl {args.wdl} \\\n"
        f"    --sample_csv {csv_path} \\\n"
        f"    --input_json_path '{json_dir}/${{SAMPLE_ID}}_hic_qc_workflow.json'"
    )
    print()


if __name__ == "__main__":
    main()
