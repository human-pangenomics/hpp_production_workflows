# T2T / findAssemblyBreakpoints

Identifies telomere-to-telomere (T2T) contigs and scaffolds and locates assembly gaps and internal telomere sequences. A contig is counted as T2T if it carries a telomere repeat array at both ends with no internal gap.

## Outputs

| File | Description |
|------|-------------|
| `t2t_ctgs`, `t2t_scfs` | Count of T2T contigs and scaffolds |
| `{assembly}.gaps.bed` | Assembly gap positions |
| `{assembly}.telomere.bed` | Telomere repeat locations |
