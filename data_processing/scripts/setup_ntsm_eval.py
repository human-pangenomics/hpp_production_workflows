#!/usr/bin/env python3
"""
Generate an input JSON for ntsm_eval_workflow and print the toil-wdl-runner command.

Collects all ntsm count files listed in a text file (one path per line) and
runs ntsmEval across all of them in a single job.

Usage:
    python setup_ntsm_eval.py count_files.txt --output_prefix my_batch
    python setup_ntsm_eval.py count_files.txt --output_prefix my_batch --partition long
"""

import argparse
import json
import os
import sys


WDL = "/private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/ntsm_eval_workflow.wdl"
TOIL_SCRIPT = "/private/nanopore/tools/hpp_production_workflows/data_processing/scripts/toil_sbatch_single_job.sh"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "file_list",
        help="Text file with one ntsm count file path per line"
    )
    parser.add_argument(
        "--output_prefix", required=True,
        help="Prefix for the output eval TSV (e.g. batch name)"
    )
    parser.add_argument(
        "--partition", default="long",
        help="Slurm partition (default: long)"
    )
    parser.add_argument(
        "--wdl", default=WDL,
        help=f"Path to ntsm_eval_workflow.wdl (default: {WDL})"
    )
    parser.add_argument(
        "--toil_script", default=TOIL_SCRIPT,
        help=f"Path to toil_sbatch_single_job.sh (default: {TOIL_SCRIPT})"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.file_list) as fh:
        files = [line.strip() for line in fh if line.strip()]

    if not files:
        sys.exit(f"No file paths found in {args.file_list}")

    cwd = os.getcwd()
    os.makedirs(os.path.join(cwd, "slurm_logs"), exist_ok=True)

    entry = {
        "ntsm_eval_wf.count_files": files,
        "ntsm_eval_wf.output_prefix": args.output_prefix,
    }

    json_path = os.path.join(cwd, f"{args.output_prefix}_ntsm_eval_inputs.json")
    with open(json_path, "w") as fh:
        json.dump(entry, fh, indent=2)

    analysis_dir = os.path.join(cwd, "analysis")
    output_file = os.path.join(cwd, f"{args.output_prefix}_ntsm_eval_outputs.json")

    print(f"\nCreated {json_path} ({len(files)} count files)\n")
    print("Run the following to submit:\n")
    print(
        f"sbatch \\\n"
        f"    --job-name={args.output_prefix}_ntsm_eval \\\n"
        f"    --partition={args.partition} \\\n"
        f"    {args.toil_script} \\\n"
        f"    --wdl {args.wdl} \\\n"
        f"    --input_json {json_path} \\\n"
        f"    --output_dir {analysis_dir} \\\n"
        f"    --output_file {output_file}"
    )
    print()


if __name__ == "__main__":
    main()
