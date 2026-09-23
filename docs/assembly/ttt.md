# TTT

## Assembly 7 — TTT (Trivial Tangle Traverser)

[**TTT**](https://github.com/marbl/TTT) (v0.1.2) resolves gaps in the raw assembly by finding traversal paths through complex regions in the assembly graph — tangles where standard graph linearization fails. This is not a polishing step; it fills sequence gaps at the graph level using read alignment evidence.

TTT targets complex repeats with at most two haplotypes per tangle. It cannot help regions with no read coverage, and does not handle rDNA tangles. The user must identify boundary nodes for each tangle manually — TTT does not auto-discover them.

!!! info "Reference"
    Antipov et al., bioRxiv 2026. "Automatic Generation of Model Sequences for Complex Regions in Assembly Graphs." [Preprint.](https://www.biorxiv.org/content/10.64898/2026.03.06.710180v1)

### Inputs

| Argument | Description |
|----------|-------------|
| `--graph` | GFA file with graph structure (or use `--verkko-output`) |
| `--alignment` | GraphAligner alignment file in GAF format |
| `--boundary-nodes` | TSV of tab-separated incoming/outgoing boundary node pairs, one pair per line. Nodes must be non-repetitive, heterozygous (for diploid tangles), and fully separate the tangle from the rest of the graph. |
| `--outdir` | Output directory |

### Running from Verkko output

```bash
./TTT.py \
  --verkko-output <verkko_output_directory> \
  --boundary-nodes boundary_nodes.tsv \
  --outdir ttt_results/
```

### Verkko ≤ v2.3.x coverage workaround

In Verkko v2.3.x, coverage of short nodes in the final graph (`assembly.homopolymer-compressed.gfa`) is unreliable. The recommended fix is to map coverage from the HiFi-only intermediate graph onto the final graph before running TTT:

```bash
# Map utig4- node IDs in the final graph to utig1- IDs in the HiFi graph
./verkko_coverage_fix/utig4_to_utig1.py <assembly_folder> > utig42utig1.gaf

# Update ONT coverage values in the final GFA using the mapped HiFi-graph coverage
./verkko_coverage_fix/utig4_coverage_updater.py \
  utig42utig1.gaf \
  <assembly_folder>/assembly.homopolymer-compressed.noseq.gfa \
  <assembly_folder>/2-processGraph/unitig-unrolled-hifi-resolved.ont-coverage.csv \
  > utig4_upt.ont-coverage.csv

# Pass corrected coverage to TTT
./TTT.py \
  --graph <assembly_folder>/assembly.homopolymer-compressed.gfa \
  --alignment <alignment.gaf> \
  --coverage utig4_upt.ont-coverage.csv \
  --boundary-nodes boundary_nodes.tsv \
  --outdir ttt_results/
```

### Outputs

| File | Description |
|------|-------------|
| `traversal.multiplicities.csv` | Node multiplicity estimates (Bandage-compatible) |
| `traversal.gaf` | Resulting traversal path in GAF format |
| `traversal.hpc.fasta` | Patch sequence (HPC-compressed for Verkko graphs; uncompressed for hifiasm) |

!!! tip "Getting non-HPC sequence from Verkko"
    Rerun Verkko with `--path traversal.gaf` — see Verkko's documentation for consensus generation from user-provided paths.

S3 prefix: `s3://human-pangenomics/submissions/4609cf7a-82af-432c-af14-9477dfbb9972--R3_verkko-v2.3.2_TTT_qc_pilot/`
