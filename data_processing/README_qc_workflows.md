# HPRC QC Workflows

## Overview

Per-sample QC for HiFi, ONT, and HiC (OmniC) data. Each workflow runs independently per file/lane and produces a summary TSV and an ntsm count file. After all samples are processed, a separate batch evaluation step runs ntsm across all count files to detect sample swaps.

---

## ntsm update (v1.2.1)

ntsm has been updated to v1.2.1. Key changes from previous versions:

- **Decoupled evaluation**: `ntsm_count` (per-sample) and `ntsm_eval` (batch) are now separate workflow steps. Per-sample QC workflows produce count files only; evaluation is submitted separately once all samples are complete. This allows flexible batch comparisons (e.g. re-running eval when new samples are added without reprocessing existing ones).
- **New reference sites**: Uses `human_sites_n10.fa` (96,287 sites from 1000 Genomes), replacing the previous sites file.
- **Coverage cap**: `ntsm_count` runs with `-m 10` (10x max coverage cap) by default to normalize across data types with different depths.
- **Clean output naming**: `ntsm_eval` symlinks count files by basename before calling `ntsmEval` so output filenames reflect sample names rather than container paths.

### Docker images

| Image | Status |
|---|---|
| `iviolich/ntsm:1.2.1` | Temporary; in use until push access to humanpangenomics org is granted |
| `humanpangenomics/ntsm:1.2.1` | Target — update WDL files once available |
| `iviolich/ont_summary_stats:1.0.1` | Temporary; fixes `calculate_summary_stats.py` bug (see below) |
| `humanpangenomics/ont_summary_stats:1.0.1` | Target |

---

## Workflows

### Per-sample QC

| Workflow | WDL | Setup script | Input format |
|---|---|---|---|
| HiFi | `hifi_qc_workflow.wdl` | `setup_hifi_qc.py` | One file path per line (BAM or fastq.gz) |
| ONT | `ont_qc_workflow.wdl` | `setup_ont_qc.py` | `reads_path<TAB>sequencing_summary_path` per line |
| HiC | `hic_qc_workflow.wdl` | `setup_hic_qc.py` | `R1_path<TAB>R2_path` per lane per line |

Each workflow:
1. Counts bases / calculates coverage
2. Runs `ntsm_count` to generate a k-mer count file per input file/lane
3. Outputs a summary TSV (coverage stats) and an ntsm count file

### Batch ntsm evaluation

After all per-sample QC jobs complete, collect all count files and run `ntsm_eval_workflow.wdl` to produce a single TSV with all pairwise comparisons. Run all data types together so cross-data-type comparisons are available.

---

## Step-by-step: running QC

### 1. Prepare a working directory

```bash
mkdir -p /path/to/qc/{hifi_blood,hifi_lcl,ont,hic,ntsm_eval}
```

Create separate subdirectories per data type. Outputs land relative to the directory from which `sbatch` is submitted (see step 4).

### 2. Build input file lists

**HiFi** — one BAM (`hifi_reads` only, not `fail_reads`) per line:
```
/path/to/HG08434/m84046_251223_003203_s1.hifi_reads.bc2077.bam
s3://human-pangenomics/submissions/.../HG08444/PacBio_HiFi/m84081_260128_015537_s3.hifi_reads.bc2004.bam
```

S3 paths are supported — Toil workers download files directly.

**ONT** — reads BAM and sequencing summary, tab-separated, one per line:
```
/path/to/HG08434_1.bam	/path/to/HG08434_1_summary.txt.gz
```

Use actual file paths, not symlinks (symlinks in upload directories may be broken).

**HiC** — R1 and R2 fastq.gz per lane, tab-separated, one lane per line:
```
/path/to/HG08434_1_LIB141063_..._R1_001.fastq.gz	/path/to/HG08434_1_LIB141063_..._R2_001.fastq.gz
```

Submit one job per lane (not per library) to preserve per-lane traceability in sample swap detection.

### 3. Generate input JSONs

Run each setup script from inside the corresponding data type directory:

```bash
source ~/toil_env/bin/activate

cd /path/to/qc/hifi_blood
python3 /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/setup_hifi_qc.py \
    files.txt --job_name hifi_blood_qc

cd /path/to/qc/ont
python3 /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/setup_ont_qc.py \
    files.txt --job_name ont_qc

cd /path/to/qc/hic
python3 /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/setup_hic_qc.py \
    files.txt --job_name hic_qc
```

Each script creates `input_jsons/`, `samples.csv`, `slurm_logs/`, and `analysis/` directories, and prints the sbatch command to run.

**Sample IDs are parsed automatically from the file path** using the regex `(HG|NA|PG)\d{5}`. The sample ID is built as `{HPRC_ID}_{filename}`. `PG` prefixes are normalized to `HG`. If no HPRC ID is found in the path, the raw filename is used and a warning is printed.

### 4. Submit jobs

**Critical**: submit each `sbatch` from inside the corresponding data type directory. The Toil script creates sample subdirectories relative to the submission working directory — submitting from the wrong directory mixes all outputs together.

```bash
cd /path/to/qc/hifi_blood && sbatch \
    --job-name=hifi_blood_qc \
    --array=[1-N]%10 \
    --partition=long \
    /private/nanopore/hprc_qc/scripts/toil_sbatch_slurm.sh \
    --wdl /private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/hifi_qc_workflow.wdl \
    --sample_csv /path/to/qc/hifi_blood/samples.csv \
    --input_json_path '/path/to/qc/hifi_blood/input_jsons/${SAMPLE_ID}_hifi_qc_workflow.json'
```

The `%10` limits concurrent array tasks to 10. Adjust as needed.

### 5. Check outputs

Outputs land in `<sample_id>/analysis/<workflow>_outputs/<uuid>/`. Each job produces:
- `<sample_id>.summary.tsv` — coverage stats
- `<sample_id>_hifi_counts.txt` (or `_ont_counts.txt`, `_hic_counts.txt`) — ntsm count file

### 6. Aggregate summary TSVs

```bash
python3 /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/aggregate_qc_summary.py \
    */analysis/*/*.summary.tsv > combined_summary.tsv
```

---

## Step-by-step: batch ntsm evaluation

Run after all per-sample QC jobs are complete.

### 1. Collect count files

```bash
QC=/path/to/qc
find $QC/hifi_blood -name "*_hifi_counts.txt" | sort >  $QC/ntsm_eval/count_files.txt
find $QC/hifi_lcl   -name "*_hifi_counts.txt" | sort >> $QC/ntsm_eval/count_files.txt
find $QC/ont        -name "*_ont_counts.txt"  | sort >> $QC/ntsm_eval/count_files.txt
find $QC/hic        -name "*_hic_counts.txt"  | sort >> $QC/ntsm_eval/count_files.txt
```

Verify the total count matches expectations before submitting.

### 2. Generate input JSON and submit

```bash
cd /path/to/qc/ntsm_eval

python3 /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/setup_ntsm_eval.py \
    count_files.txt --output_prefix <batch_name>

source ~/toil_env/bin/activate
cd /path/to/qc/ntsm_eval && sbatch \
    --job-name=<batch_name>_ntsm_eval \
    --partition=long \
    /private/nanopore/tools/hpp_production_workflows/data_processing/scripts/toil_sbatch_single_job.sh \
    --wdl /private/nanopore/tools/hpp_production_workflows/data_processing/wdl/workflows/ntsm_eval_workflow.wdl \
    --input_json /path/to/qc/ntsm_eval/<batch_name>_ntsm_eval_inputs.json \
    --output_dir /path/to/qc/ntsm_eval/analysis \
    --output_file /path/to/qc/ntsm_eval/<batch_name>_ntsm_eval_outputs.json
```

### 3. Interpret results

Output: `analysis/<uuid>/<batch_name>_ntsm_eval.tsv`

Key columns:
- `same` — 1 if samples match, 0 if swap detected
- `relate` — relatedness score (>0.98 expected for same sample across data types)
- `ibs0` — sites where both samples are homozygous for different alleles (should be ~0 for same sample)

---

## Notes

- **Toil cleanup**: Do not run `toil clean` on a job store while outputs are still needed.
- **ONT summary files**: Use actual file paths, not symlinks. Symlinks in upload directories may point to pre-basecalling paths that no longer exist.
- **HiFi methylation**: The HiFi workflow includes a methylation check by default. Pass `--no_methylation` to `setup_hifi_qc.py` to skip it.
- **Kinnex**: Not yet supported. Will require a custom ntsm sites file. TBD.
- **Docker tags**: Update WDL files from `iviolich/ntsm:1.2.1` → `humanpangenomics/ntsm:1.2.1` and `iviolich/ont_summary_stats:1.0.1` → `humanpangenomics/ont_summary_stats:1.0.1` once push access is available.
- **ONT `calculate_summary_stats.py` bug**: Fixed in `iviolich/ont_summary_stats:1.0.1`. Original script crashed with `IndexError` when the fail-reads file contained no reads. Both empty and single-read cases are now handled.

---

## Troubleshooting

**Jobs fail with `coordinationDir (/data/tmp) does not exist`**  
A cluster node has `/data/tmp` unavailable. Exclude it and resubmit:
```bash
sbatch --exclude=<node_name> ...
```

**Jobs fail with `job store already exists`**  
A previous run left a jobstore in the sample directory. If the failure was infrastructure-related, resume with:
```bash
sbatch ... --toil_args '--restart'
```
To start fresh: `rm -rf <sample_id>/<workflow>_jobstore`

**ONT jobs fail in `calc_ont_summary_stats` but count files are present**  
The ntsm count task ran before the stats task failed. The count file is in the jobstore at `<sample_id>/ont_qc_workflow_jobstore/files/for-job/.../`. Fix the underlying issue and rerun with `--restart`, or rescue the file manually.
