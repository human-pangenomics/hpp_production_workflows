# Panpatch

## Assembly 8 — Panpatch

Panpatch is a self-reference assembly patching step applied to a pilot cohort of TTT assemblies. Each sample's own HiFi and ONT reads are used to fill gaps and improve continuity at regions that remain unresolved after TTT.

Pilot cohort (**panpatch pilot-v1-ttt**): 18 samples. Patched assemblies:

```
s3://human-pangenomics/submissions/8b47ec5e-9653-11f1-8a5c-e76d6dc4bd3d--panpatch-pilot-v1/pilot-v1-ttt/patches-v2/{sid}.hap{1,2}.fa.gz
```

After patching, apply the full QC pipeline to the patched FASTAs. Two steps differ from a standard QC run:

- **CenSat** must be re-run — the existing BED from the TTT assembly is invalid after patching (see [QC 2 · CenSat](../qc/censat.md))
- **Flagger** has a different read policy — DC FASTQs supersede raw HiFi BAMs where available (see [QC 7 · Flagger](../qc/flagger.md))
