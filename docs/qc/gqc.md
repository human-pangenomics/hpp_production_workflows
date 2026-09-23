# GQC

Validates the assembly against the HG002 v1.1 benchmark using GQC. Download GQC and the benchmark from the [GQC repository](https://github.com/marbl/GQC) first.

### 1 — Point at the config and benchmark

```bash
CONFIGFILE=<path to GQC>/GQC/benchconfig.txt

ln -s <path to v1.1>/v1.1.fasta.gz     || true
ln -s <path to v1.1>/v1.1.fasta.gz.gzi || true
ln -s <path to v1.1>/v1.1.fasta.gz.fai || true
```

The `|| true` keeps the block from aborting if symlinks already exist from a previous run.

### 2 — Locate the assembly

```bash
asm="assembly.fasta"
if [ ! -e "$asm" ]; then
  asm="assembly.fasta.gz"
fi
if [ ! -e "$asm" ]; then
  echo "Error: cannot find either assembly.fasta or assembly.fasta.gz"
  exit 1
fi
```

### 3 — Run GQC

```bash
GQC -a minimap2 -c $CONFIGFILE -r v1.1.fasta.gz -q $asm \
  -p gqc -A verkko_hg002 -B v1.1
```

`-r` is the reference benchmark and `-q` the assembly under test. `-A` and `-B` are the labels applied to each in the report.
