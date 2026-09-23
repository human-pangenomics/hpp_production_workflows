# Verkko Assembly

## Assembly 4 — Assembly with verkko

Requires verkko v2.3 or later.

```bash
verkko --slurm -d verkko-hi-c \
  --ovb-run 8 32 32 \
  --screen-human-contaminants \
  --hifi $HIFI r10-hifiasm-correct.ec.ont.fq.gz \
  --nano $ONT $ONT_R10 \
  --hic1 $HIC1 \
  --hic2 $HIC2 \
  --unitig-abundance 8
```

The corrected ONT UL is passed to `--hifi` alongside the raw HiFi, because correction has brought it to high-quality-sequence accuracy. The raw ONT UL goes to `--nano` as resolving data, so the same underlying reads appear on both sides in different forms. This is intentional, not a copy-paste error.

!!! note "`--unitig-abundance 8`"
    Recommended for high-coverage datasets. Combined ONT UL and HiFi coverage here is 180×.

!!! note "`--slurm`"
    Auto-requests memory and time per job. Drop it to run on a single node, where verkko will auto-detect available CPUs and memory.

Expect around 41 hours walltime on a cluster. Verkko checkpoints internally, so re-running the same command in the same `-d` directory resumes rather than restarting.

---

## Assembly 5 — Polishing

No polishing step is required. Assembly QV is estimated at 56 using yak with k=31 mers from Illumina data, or Q48 as measured by GQC. Polishing at this accuracy is not recommended and is unlikely to improve the result.

---

## Checkpoints

Worth confirming before moving to the next stage:

| After | Check |
|-------|-------|
| 3.3 | ONT input is smaller than the raw input — the 20 kbp filter should remove reads |
| 3.4 | `--hom-cov` lands near your expected total coverage; a wildly off value usually means a missing input file |
| 3.7 | Both ID files are non-empty and the HiFi/ONT ratio is plausible |
| 3.8 | `r10-hifiasm-correct.ec.ont.fq.gz` opens cleanly before deleting anything |
| 4 | Assembly size is near 6.2 Gbp across both haplotypes |
