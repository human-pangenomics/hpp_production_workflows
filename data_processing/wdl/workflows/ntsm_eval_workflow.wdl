version 1.0

import "../tasks/ntsm.wdl" as ntsm_tasks

workflow ntsm_eval_wf {
    input {
        Array[File] count_files
        String output_prefix
    }

    parameter_meta {
        count_files: "All ntsm count files to evaluate together (from ntsm_count task across all samples and read types)."
        output_prefix: "Prefix for the output eval TSV (e.g. batch name or project name)."
    }

    call ntsm_tasks.ntsm_eval {
        input:
            count_files   = count_files,
            output_prefix = output_prefix
    }

    call summarize_ntsm_eval {
        input:
            ntsm_eval_tsv = ntsm_eval.ntsm_eval_out,
            output_prefix = output_prefix
    }

    output {
        File ntsm_eval_out     = ntsm_eval.ntsm_eval_out
        File ntsm_eval_summary = summarize_ntsm_eval.summary_tsv
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Runs ntsmEval across all ntsm count files from a batch of QC runs to detect sample swaps."
    }
}


task summarize_ntsm_eval {
    input {
        File   ntsm_eval_tsv
        String output_prefix

        Int memSizeGB   = 4
        Int threadCount = 1
        Int disk_size   = 16
        Int preempts    = 2
    }

    command <<<
set -euo pipefail

cat << 'PYEOF' > /tmp/summarize_ntsm.py
import sys, csv, re
from collections import defaultdict

tsv_path    = sys.argv[1]
output_prefix = sys.argv[2]

HPRC_ID_RE = re.compile(r'((?:HG|NA|PG)\d{5})')

def get_sample(name):
    m = HPRC_ID_RE.search(name)
    return m.group(1) if m else None

def get_type(name):
    if 'hifi'   in name: return 'hifi'
    if 'ont'    in name: return 'ont'
    if 'hic'    in name: return 'hic'
    if 'kinnex' in name: return 'kinnex'
    return 'other'

files  = defaultdict(lambda: defaultdict(set))
within = defaultdict(list)
cross  = defaultdict(list)
all_samples = set()

with open(tsv_path) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        s1 = get_sample(row['sample1'])
        s2 = get_sample(row['sample2'])
        score = float(row['score'])
        if s1:
            files[s1][get_type(row['sample1'])].add(row['sample1'])
            all_samples.add(s1)
        if s2:
            files[s2][get_type(row['sample2'])].add(row['sample2'])
            all_samples.add(s2)
        if s1 is None or s2 is None:
            continue
        if s1 == s2:
            within[s1].append(score)
        else:
            cross[s1].append(score)
            cross[s2].append(score)

samples = sorted(all_samples)
all_types = ['hifi', 'ont', 'hic', 'kinnex', 'other']
type_cols = [t for t in all_types if any(files[s][t] for s in samples)]

total_files = sum(sum(len(files[s][t]) for t in type_cols) for s in samples)

cols = (
    ['sample']
    + [f'n_{t}' for t in type_cols]
    + ['n_files',
       'n_within', 'exp_within', 'min_within(exp>0)', 'mean_within(exp<0.25)', 'max_within(exp<0.25)',
       'n_cross',  'exp_cross',  'min_cross(exp>0.5)', 'mean_cross(exp>0.5)',  'max_cross(exp<3.0)',
       'flag']
)

rows = []
for s in samples:
    ws = within[s]
    cs = cross[s]
    if not ws or not cs:
        print(f"Warning: no data for {s}", file=sys.stderr)
        continue
    min_w  = min(ws)
    mean_w = sum(ws) / len(ws)
    max_w  = max(ws)
    min_c  = min(cs)
    mean_c = sum(cs) / len(cs)
    max_c  = max(cs)
    flag   = 'PASS' if min_w > 0 and max_w < 0.5 and min_c > 0.5 else 'FAIL'
    n_by_type = [len(files[s][t]) for t in type_cols]
    n = sum(n_by_type)
    exp_w = n * (n - 1) // 2
    exp_c = n * (total_files - n)
    rows.append(
        [s]
        + n_by_type
        + [n,
           len(ws), exp_w, f'{min_w:.3f}', f'{mean_w:.3f}', f'{max_w:.3f}',
           len(cs), exp_c, f'{min_c:.3f}', f'{mean_c:.3f}', f'{max_c:.3f}',
           flag]
    )

out_file = output_prefix + '_ntsm_sample_summary.tsv'
with open(out_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(cols)
    writer.writerows(rows)

n_fail = sum(1 for r in rows if r[-1] == 'FAIL')
print(f"Samples: {len(rows)} | PASS: {len(rows)-n_fail} | FAIL: {n_fail}")
print(f"Written to {out_file}")
PYEOF

python3 /tmp/summarize_ntsm.py ~{ntsm_eval_tsv} ~{output_prefix}
    >>>

    output {
        File summary_tsv = "~{output_prefix}_ntsm_sample_summary.tsv"
    }

    runtime {
        memory:      memSizeGB + " GB"
        cpu:         threadCount
        disks:       "local-disk " + disk_size + " SSD"
        docker:      "python:3.11-slim"
        preemptible: preempts
    }
}
