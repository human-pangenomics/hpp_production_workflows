#!/usr/bin/env python3
"""
Generate input JSONs and samples CSV for ont_qc_workflow, then print the sbatch command.

Usage:
    python setup_ont_qc.py files.txt
    python setup_ont_qc.py files.txt --job_name my_batch --partition long

Input file format (tab-separated, one sample per line):
    /path/to/sample1.fastq.gz    /path/to/sample1_sequencing_summary.txt
    /path/to/sample2.bam         /path/to/sample2_sequencing_summary.txt
"""

import argparse
import json
import os
import re
import csv
import sys


WDL = "/private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/ont_qc_workflow.wdl"
TOIL_SCRIPT = "/private/nanopore/hprc_qc/scripts/toil_sbatch_slurm.sh"

HPRC_ID_RE = re.compile(r'((?:HG|NA|PG)\d{5})')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "file_list",
        help="Tab-separated text file with one sample per line: reads_path<TAB>sequencing_summary_path"
    )
    parser.add_argument(
        "--job_name", default="ont_qc",
        help="Slurm job name (default: ont_qc)"
    )
    parser.add_argument(
        "--partition", default="long",
        help="Slurm partition (default: long)"
    )
    parser.add_argument(
        "--wdl", default=WDL,
        help=f"Path to ont_qc_workflow.wdl (default: {WDL})"
    )
    return parser.parse_args()


def strip_ext(path):
    name = os.path.basename(path)
    for ext in (".bam", ".fastq.gz", ".fastq", ".fq.gz", ".fq"):
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
            sys.exit(f"Line {i} must have exactly two tab-separated fields (reads<TAB>summary):\n  {line}")
        reads_path, summary_path = parts
        samples.append((reads_path.strip(), summary_path.strip()))

    cwd = os.getcwd()
    json_dir = os.path.join(cwd, "input_jsons")
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(os.path.join(cwd, "slurm_logs"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "analysis"), exist_ok=True)

    csv_rows = []
    for reads_path, summary_path in samples:
        sid = build_sample_id(reads_path)
        entry = {
            "ont_qc_wf.ont_reads": reads_path,
            "ont_qc_wf.ont_sequencing_summary_files": summary_path,
            "ont_qc_wf.file_id": sid,
        }
        json_path = os.path.join(json_dir, f"{sid}_ont_qc_workflow.json")
        with open(json_path, "w") as fh:
            json.dump(entry, fh, indent=2)
        csv_rows.append({"sample_id": sid, "ont_reads": reads_path, "sequencing_summary": summary_path})

    csv_path = os.path.join(cwd, "samples.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "ont_reads", "sequencing_summary"])
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
        f"    --input_json_path '{json_dir}/${{SAMPLE_ID}}_ont_qc_workflow.json'"
    )
    print()


if __name__ == "__main__":
    main()
