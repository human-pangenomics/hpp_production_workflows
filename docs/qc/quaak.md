# Quaak

Unique k-mer completeness assessment benchmarked against CHM13v2. Builds synteny blocks between the assembly and reference to measure how much of the unique k-mer space is covered. Regions where expected unique k-mers are absent or misplaced flag potential assembly errors or large SVs.

## Outputs

| File | Metric |
|------|--------|
| `*.kma` | k-mer match accuracy — fraction of CHM13v2 unique k-mers found in the assembly |
| `*.ksa_pct` | k-mer synteny accuracy — fraction found in syntenic positions |
