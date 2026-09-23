# Compleasm

BUSCO-based gene completeness assessment using the `primates_odb10` database. Reports the fraction of single-copy orthologs that are complete (single copy or duplicated), fragmented, or missing. Run per haplotype FASTA.

```bash
compleasm run \
  -a {haplotype.fasta} \
  -o {outdir} \
  -l primates_odb10 \
  -t {threads}
```

## Outputs

| File | Description |
|------|-------------|
| `summary.txt` | S/D/F/M gene counts and percentages |
