# Data Processing

Per-sample QC for HiFi, ONT, and Hi-C data. Each workflow runs independently per file/lane and produces a summary TSV and an ntsm count file. After all samples are processed, a separate batch evaluation step runs ntsm across all count files to detect sample swaps.

## Workflows

| Workflow | WDL | Input format |
|----------|-----|--------------|
| HiFi QC | `hifi_qc_workflow.wdl` | One file path per line (BAM or fastq.gz) |
| ONT QC | `ont_qc_workflow.wdl` | `reads_path<TAB>sequencing_summary_path` per line |
| Hi-C QC | `hic_qc_workflow.wdl` | `R1_path<TAB>R2_path` per lane per line |

Each workflow:

1. Counts bases / calculates coverage
2. Runs `ntsm_count` to generate a k-mer count file per input file/lane
3. Outputs a summary TSV (coverage stats) and an ntsm count file

## Batch ntsm evaluation

After all per-sample QC jobs complete, collect all count files and run `ntsm_eval_workflow.wdl` to produce a single TSV with all pairwise comparisons. Run all data types together so cross-data-type comparisons are available.

## ntsm v1.2.1 changes

- **Decoupled evaluation**: `ntsm_count` (per-sample) and `ntsm_eval` (batch) are now separate workflow steps
- **New reference sites**: Uses `human_sites_n10.fa` (96,287 sites from 1000 Genomes)
- **Coverage cap**: `ntsm_count` runs with `-m 10` (10× max coverage cap) to normalize across data types
- **Clean output naming**: `ntsm_eval` symlinks count files by basename so output filenames reflect sample names

WDL source: [`data_processing/wdl/`](https://github.com/human-pangenomics/hpp_production_workflows/tree/v3/data_processing/wdl)
