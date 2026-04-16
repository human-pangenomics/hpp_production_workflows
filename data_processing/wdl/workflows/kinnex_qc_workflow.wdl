version 1.0

import "../tasks/ntsm.wdl" as ntsm_check

workflow kinnex_qc_wf {
    input {
        File   flnc_bam
        File   pigeon_summary
        String sample_id
    }

    parameter_meta {
        flnc_bam:       "FLNC BAM produced by SMRTLink IsoSeq pipeline (refine step)."
        pigeon_summary: "pigeon_filtered.summary.txt produced by SMRTLink IsoSeq pipeline."
        sample_id:      "Sample ID used for output file naming."
    }

    meta {
        author: "Ivo Violich"
        email:  "iviolich@ucsc.edu"
    }

    ## ntsm k-mer count from FLNC reads for cross-sample swap detection
    call ntsm_check.ntsm_count as ntsm_wf {
        input:
            input_reads = [flnc_bam],
            sample_id   = sample_id,
            read_type   = "kinnex"
    }

    ## parse pigeon classification summary into a single-row TSV
    call summarize_kinnex_qc {
        input:
            pigeon_summary = pigeon_summary,
            file_name      = sample_id
    }

    output {
        ## per-sample QC summary
        File kinnex_qc_summary  = summarize_kinnex_qc.summary_file

        ## ntsm count file (batch eval run separately across all samples)
        File kinnex_ntsm_counts = ntsm_wf.ntsm_counts
    }
}


task summarize_kinnex_qc {
    input {
        File   pigeon_summary
        String file_name

        Int memSizeGB   = 4
        Int threadCount = 1
        Int disk_size   = 16
        Int preempts    = 2
    }

    command <<<
set -euo pipefail

## Write parser to a temp file (avoids WDL/bash interpolation inside heredoc)
cat << 'PYEOF' > /tmp/parse_pigeon.py
import sys
import re

pigeon_file = sys.argv[1]
file_name   = sys.argv[2]

def parse_count_pct(line):
    m = re.search(r':\s+(\d+)\s+\(([0-9.]+)%\)', line)
    if m:
        return m.group(1), m.group(2)
    m = re.search(r':\s+(\d+)', line)
    if m:
        return m.group(1), 'NA'
    return 'NA', 'NA'

with open(pigeon_file) as f:
    content = f.read()

data = {}

for line in content.splitlines():
    s = line.strip()
    if s.startswith('Input'):
        data['flnc_input_reads'] = s.split(':')[1].strip()
    elif s.startswith('Passed'):
        m = re.search(r':\s+(\d+)\s+\(([0-9.]+)%\)', s)
        if m:
            data['passed_reads'] = m.group(1)
            data['passed_pct']   = m.group(2)
    elif s.startswith('Unique genes'):
        data['unique_genes'] = s.split(':')[1].strip()
    elif s.startswith('Unique transcripts'):
        data['unique_transcripts'] = s.split(':')[1].strip()

## by-isoform classification counts (unique isoforms per category)
if 'Classifications, by isoform' in content and 'Classifications, by read' in content:
    by_isoform = content.split('Classifications, by isoform')[1].split('Classifications, by read')[0]
    for line in by_isoform.splitlines():
        s = line.strip()
        if s.startswith('Full splice match'):
            data['fsm_isoforms'], data['fsm_isoform_pct'] = parse_count_pct(s)
        elif s.startswith('Incomplete splice match'):
            data['ism_isoforms'], data['ism_isoform_pct'] = parse_count_pct(s)
        elif s.startswith('Novel in catalog'):
            data['nic_isoforms'], data['nic_isoform_pct'] = parse_count_pct(s)
        elif s.startswith('Novel not in catalog'):
            data['nnc_isoforms'], data['nnc_isoform_pct'] = parse_count_pct(s)

## by-read FSM % (fraction of reads mapping to known isoforms)
if 'Classifications, by read' in content:
    by_read = content.split('Classifications, by read')[1].split('Junctions')[0]
    for line in by_read.splitlines():
        s = line.strip()
        if s.startswith('Full splice match'):
            _, data['fsm_read_pct'] = parse_count_pct(s)
            break

## filter reasons
if 'Filter reasons' in content:
    filter_section = content.split('Filter reasons')[1]
    for line in filter_section.splitlines():
        s = line.strip()
        if s.startswith('Intrapriming'):
            _, data['intrapriming_pct'] = parse_count_pct(s)
        elif s.startswith('RT switching'):
            _, data['rt_switching_pct'] = parse_count_pct(s)

cols = [
    'file_name',
    'flnc_input_reads', 'passed_reads', 'passed_pct',
    'unique_genes', 'unique_transcripts',
    'fsm_isoforms', 'fsm_isoform_pct',
    'ism_isoforms', 'ism_isoform_pct',
    'nic_isoforms', 'nic_isoform_pct',
    'nnc_isoforms', 'nnc_isoform_pct',
    'fsm_read_pct',
    'intrapriming_pct', 'rt_switching_pct',
]

row = [file_name] + [data.get(c, 'NA') for c in cols[1:]]

out_file = file_name + '.summary.tsv'
with open(out_file, 'w') as out:
    out.write('\t'.join(cols) + '\n')
    out.write('\t'.join(row) + '\n')

print(f"Generated: {out_file}")
PYEOF

python3 /tmp/parse_pigeon.py ~{pigeon_summary} ~{file_name}
    >>>

    output {
        File summary_file = "~{file_name}.summary.tsv"
    }

    runtime {
        memory:     memSizeGB + " GB"
        cpu:        threadCount
        disks:      "local-disk " + disk_size + " SSD"
        docker:     "python:3.11-slim"
        preemptible: preempts
    }
}
