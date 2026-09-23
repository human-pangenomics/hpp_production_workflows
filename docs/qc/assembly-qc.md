# assembly_qc / Yak QV

k-mer quality value estimation (k=31) using Hi-C reads as the k-mer database. Yak counts k-mers in the Hi-C reads and compares them against the assembly to estimate per-base error rate and overall QV. Run via the `standard_qc_nontrio` WDL workflow.

!!! warning "Memory"
    ~500 GB RAM required per sample. Schedule accordingly on high-memory nodes.

## Outputs

| File | Description |
|------|-------------|
| `{assembly}.yak.qv` | Per-haplotype QV score |
