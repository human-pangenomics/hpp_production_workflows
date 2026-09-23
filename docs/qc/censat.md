# CenSat

Centromere and satellite sequence annotation using RepeatMasker. Produces a per-haplotype BED of all satellite repeat classes and a merged diploid BED.

The diploid BED (`cenSat.dip.bed`) is a required input to Flagger — it masks centromeric regions during read coverage modeling, preventing false HMM state calls in regions where coverage depth is biologically extreme.

Run as WDL (`cenSat.wdl`) via Toil on SLURM.

!!! warning "After Panpatch"
    The existing CenSat BED from the TTT assembly is invalid after patching — the sequence has changed. Re-run CenSat from the patched FASTA before running Flagger.

## Outputs

| File | Description |
|------|-------------|
| `{sample}.cenSat.bed` | Per-haplotype satellite annotation |
| `{sample}.cenSat.dip.bed` | Diploid merged BED — required input to Flagger |
