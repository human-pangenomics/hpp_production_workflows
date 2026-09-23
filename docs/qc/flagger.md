# Flagger

HMM-based read coverage classifier. Maps both HiFi and ONT reads to the assembly and models coverage in 16 kb windows, assigning each window one of four states:

| State | Meaning |
|-------|---------|
| **Err** | Error / no coverage |
| **Dup** | Duplicated |
| **Col** | Collapsed |
| **Hap** | Haploid / expected |

CenSat masking prevents centromeric high-coverage windows from pulling the model.

Run as WDL (`hmm_flagger_end_to_end_with_mapping.wdl`) via Toil on SLURM, one array task per sample per read type (HiFi + ONT = 2 tasks per sample).

!!! warning "Depends on CenSat"
    CenSat must complete before Flagger can start.

## Outputs

| File | Description |
|------|-------------|
| `*.hifi_flagger_prediction.bed` | HiFi-based window classifications |
| `*.ont_flagger_prediction.bed` | ONT-based window classifications |
| `*.bam`, `*.bam.bai` | Read alignments used for coverage modeling — input to NucFlag |

!!! note "BAM archival"
    Alignment BAMs are ~100–200 GB per sample. Archive to S3 after NucFlag completes.

## After Panpatch — read policy

Where DeepConsensus (DC) FASTQs are available for a HiFi movie, the DC FASTQ supersedes the raw HiFi BAM for that movie. DC data has better basecalling and more passing reads. Match DC FASTQs to HiFi movies by run ID (`m\d+_\d+_\d+`).

Also use `USE_SHARED_JOBSTORE=1` to avoid node-local scratch contention when running many samples.
