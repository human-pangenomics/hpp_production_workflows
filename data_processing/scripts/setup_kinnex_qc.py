#!/usr/bin/env python3
"""
Generate input JSONs and samples CSV for kinnex_qc_workflow, then print the sbatch command.

Usage:
    python setup_kinnex_qc.py files.txt
    python setup_kinnex_qc.py files.txt --job_name kinnex_qc --partition long

Input file format (tab-separated, one sample per line):
    /path/to/HG08434.flnc.bam\t/path/to/HG08434.pigeon_filtered.summary.txt

S3 paths are supported — Toil workers download files directly.
"""

import argparse
import json
import os
import re
import csv
import sys


WDL = "/private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/kinnex_qc_workflow.wdl"
TOIL_SCRIPT = "/private/nanopore/hprc_qc/scripts/toil_sbatch_slurm.sh"

HPRC_ID_RE = re.compile(r'((?:HG|NA|PG)\d{5})')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "file_list",
        help="Tab-separated file: flnc_bam_path<TAB>pigeon_summary_path, one sample per line"
    )
    parser.add_argument(
        "--job_name", default="kinnex_qc",
        help="Slurm job name (default: kinnex_qc)"
    )
    parser.add_argument(
        "--partition", default="long",
        help="Slurm partition (default: long)"
    )
    parser.add_argument(
        "--wdl", default=WDL,
        help=f"Path to kinnex_qc_workflow.wdl (default: {WDL})"
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
        lines = [line.strip() for line in fh if line.strip()]

    if not lines:
        sys.exit(f"No entries found in {args.file_list}")

    entries = []
    for i, line in enumerate(lines, 1):
        parts = line.split('\t')
        if len(parts) != 2:
            sys.exit(
                f"Line {i}: expected two tab-separated fields (flnc_bam<TAB>pigeon_summary), got: {line!r}"
            )
        entries.append((parts[0].strip(), parts[1].strip()))

    cwd = os.getcwd()
    json_dir = os.path.join(cwd, "input_jsons")
    os.makedirs(json_dir, exist_ok=True)
    os.makedirs(os.path.join(cwd, "slurm_logs"), exist_ok=True)
    os.makedirs(os.path.join(cwd, "analysis"), exist_ok=True)

    samples = []
    for flnc_bam, pigeon_summary in entries:
        sid = build_sample_id(flnc_bam)
        entry = {
            "kinnex_qc_wf.flnc_bam":       flnc_bam,
            "kinnex_qc_wf.pigeon_summary":  pigeon_summary,
            "kinnex_qc_wf.sample_id":       sid,
        }
        json_path = os.path.join(json_dir, f"{sid}_kinnex_qc_workflow.json")
        with open(json_path, "w") as fh:
            json.dump(entry, fh, indent=2)
        samples.append({"sample_id": sid, "flnc_bam": flnc_bam, "pigeon_summary": pigeon_summary})

    csv_path = os.path.join(cwd, "samples.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample_id", "flnc_bam", "pigeon_summary"])
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
        f"    --input_json_path '{json_dir}/${{SAMPLE_ID}}_kinnex_qc_workflow.json'"
    )
    print()


if __name__ == "__main__":
    main()
