# QC Pipeline

The QC pipeline is a standalone suite that can be applied to any assembly. Run it after verkko (Assembly 4), after TTT (Assembly 7), after Panpatch (Assembly 8), or any combination. For R3 we ran it after verkko and again after TTT.

**CenSat must complete before Flagger can start.** All other tools run in parallel with CenSat.

```mermaid
flowchart TD
  asm{{"assembly FASTA"}}
  hifi{{"HiFi reads"}}
  ont{{"ONT reads"}}
  hic{{"Hi-C reads"}}
  bench{{"HG002 v1.1 benchmark"}}

  asm & bench --> gqc["GQC"]
  asm --> q1["CenSat"]
  asm --> q2["Compleasm"]
  asm --> q3["T2T"]
  asm --> q4["Quaak"]
  asm --> q5["assembly_qc"]
  hic --> q5

  gqc --> gout(["benchmark report"])
  q1 --> cbed(["cenSat.dip.bed"])
  cbed & hifi & ont --> q6["Flagger"]

  q2 --> r2(["BUSCO summary"])
  q3 --> r3(["T2T counts + gaps BED"])
  q4 --> r4(["kma / ksa_pct"])
  q5 --> r5(["QV score"])

  q6 --> bams(["prediction BEDs + BAMs"])
  bams --> q7["NucFlag"]
  q7 --> nf(["nucflag.bed"])
```

## Tool summary

| Step | Tool | Input | Key output | Depends on |
|------|------|-------|-----------|------------|
| QC 1 | GQC | assembly + HG002 v1.1 benchmark | benchmark report | — |
| QC 2 | CenSat | assembly FASTA | `cenSat.dip.bed` | — |
| QC 3 | Compleasm | assembly FASTA | BUSCO summary | — |
| QC 4 | T2T | assembly FASTA | T2T counts, gaps BED | — |
| QC 5 | Quaak | assembly FASTA | kma / ksa_pct | — |
| QC 6 | assembly_qc / Yak QV | assembly + Hi-C reads | QV score | — |
| QC 7 | Flagger | assembly + HiFi + ONT + CenSat BED | prediction BEDs + BAMs | CenSat |
| QC 8 | NucFlag | Flagger BAMs | nucflag.bed | Flagger |
