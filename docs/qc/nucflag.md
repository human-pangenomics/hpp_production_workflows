# NucFlag

Per-contig misassembly detection using the read alignments produced by Flagger. NucFlag scans alignment depth and split-read signals to flag structural errors within contigs that Flagger's window-level HMM cannot resolve.

!!! warning "Depends on Flagger"
    Requires Flagger BAMs as input.

## Error classes

| Label | Meaning |
|-------|---------|
| `MISJOIN` | Coverage drop within a contig — two distinct sequences joined incorrectly |
| `COLLAPSE_OTHER` | Collapsed repeat from a different genomic context |
| `COLLAPSE_VAR` | Collapsed heterozygous variant — two haplotypes merged into one sequence |

## Outputs

| File | Description |
|------|-------------|
| `*.nucflag.bed` | Annotated misassembly loci with error class labels |
