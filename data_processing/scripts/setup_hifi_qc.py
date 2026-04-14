#!/usr/bin/env python3
"""
Generate input JSONs and samples CSV for hifi_qc_workflow, then print the sbatch command.

Usage:
    python setup_hifi_qc.py files.txt
    python setup_hifi_qc.py files.txt --no_methylation
    python setup_hifi_qc.py files.txt --job_name my_batch --partition long

Input file format (one path per line):
    /path/to/sample1.bam
    /path/to/sample2.bam
    /path/to/sample3.fastq.gz
"""

import argparse
import json
import os
import re
import csv
import sys


WDL = "/private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/hifi_qc_workflow.wdl"
TOIL_SCRIPT = "/private/nanopore/hprc_qc/scripts/toil_sbatch_slurm.sh"

HPRC_ID_RE = re.compile(r'((?:HG|NA|PG)\d{5})')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "file_list",
        help="Text file with one input file path per line"
    )
    parser.add_argument(
        "--no_methylation", action="store_true",
        help="Skip methylation check"
    )
    parser.add_argument(
        "--job_name", default="hifi_qc",
        help="Slurm job name (default: hifi_qc)"
    )
    parser.add_argument(
        "--partition", default="long",
        help="Slurm partition (default: long)"
    )
    parser.add_argument(
        "--wdl", default=WDL,
        help=f"Path to hifi_qc_workflow.wdl (default: {WDL})"
    )
    return parser.parse_args()


def strip_ext(path):
    name = os.path.basename(path)
    for ext in (".bam", ".fastq.gz", ".fastq", ".fq.gz", ".fq"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def normalize_sample_id(hprc_id):
    """Normalize PG##### -> HG##### (same individuals, different prefix convention)."""
    if hprc_id.startswith("PG"):
        return "HG" + hprc_id[2:]
    return hprc_id


def build_sample_id(path):
    """Prepend HPRC sample ID (HG/NA/PG#####) to file basename if found in path."""
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
        files = [line.strip() for line in fh if line.strip()]

    if not files:
        sys.exit(f"No file paths found in {args.file_list}")

    cwd = os.getcwd()
    json_dir = os.path.join(cwd, "input_jsons")
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(os.path.join(cwd, "slurm_logs"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "analysis"), exist_ok=True)

    samples = []
    for f in files:
        sid = build_sample_id(f)
        entry = {
            "hifi_qc_wf.hifi_reads": f,
            "hifi_qc_wf.sample_id": sid,
            "hifi_qc_wf.perform_methylation_check": not args.no_methylation,
        }
        json_path = os.path.join(json_dir, f"{sid}_hifi_qc_workflow.json")
        with open(json_path, "w") as fh:
            json.dump(entry, fh, indent=2)
        samples.append({"sample_id": sid, "hifi_reads": f})

    csv_path = os.path.join(cwd, "samples.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "hifi_reads"])
        writer.writeheader()
        writer.writerows(samples)

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
        f"    --input_json_path '{json_dir}/${{SAMPLE_ID}}_hifi_qc_workflow.json'"
    )
    print()


if __name__ == "__main__":
    main()
