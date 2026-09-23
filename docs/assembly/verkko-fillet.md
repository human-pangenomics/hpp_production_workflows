# verkko-fillet

## Assembly 6 — verkko-fillet

[**verkko-fillet**](https://github.com/marbl/verkko-fillet) is a Python toolkit for cleaning, fixing, and gap-filling assemblies produced by Verkko. It is used here primarily for its **internal telomere detection** capability: identifying telomere repeat arrays that fall within a contig rather than at its ends. These internal telomeres mark assembly tangles where distinct chromosomal arms have been incorrectly joined, and the resulting annotations inform the boundary node selection for TTT (Stage 7).

The toolkit is designed for **interactive Jupyter notebook workflows** rather than a fixed command-line pipeline. Key notebooks: *Verkko QC* (quality metrics, T2T status, chromosome assignment) and *Recovering T2T contigs* (internal telomere detection, contig trimming, fixed GAF path output for consensus rebuilding).

!!! info "Reference"
    Juhi Kim et al., bioRxiv 2025. "Finishing a complete giraffe genome from telomere to telomere with Verkko-Fillet."

### Installation

```bash
pip install verkkofillet
```

External dependencies required on `$PATH`: `mashmap`, `samtools` / `bgzip`, `seqtk`.

### Key outputs

| Output | Description |
|--------|-------------|
| Internal telomere annotations | Per-contig summary of internal telomere positions — used to identify TTT boundary node candidates |
| Fixed GAF path file | Corrected traversal path for consensus rebuilding with Verkko |
| QC metrics | N50, T2T status, chromosome assignment, completeness |

Can run in parallel with the QC pipeline; both take the raw verkko assembly as input and have no dependency on each other.
