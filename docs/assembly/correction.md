# Download & Correction

## Assembly 1 — Download the input data

### HiFi Revio · ~60× · BAM

PacBio Revio HiFi reads, aligned to HG002 v1.0. Used as the error-correction guide in Stage 3 and passed directly to verkko in Stage 4.

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/hifi_revio_pbmay24/hg002v1.0.1_hifi_revio_pbmay24.bam
```

### ONT R10 UL · ~60× · BAM

Oxford Nanopore R10 ultra-long reads basecalled with Dorado, aligned to HG002 v1.0. These are the reads that will be corrected in Stage 3 using the HiFi data.

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.0/mapping/ont_r10_ul_dorado/hg002v1.0_ont_r10_ul_dorado.bam
```

### ONT R9 UL · optional · FASTQ.GZ

Oxford Nanopore R9 ultra-long reads (guppy 6.3.7 rebasecalling). Not corrected; enters verkko at Stage 4 as additional resolving data. Only the three files below were used in the reference assembly. If you skip R9, leave the `ONT` variable empty in Stage 2 — the verkko command is unchanged.

```bash
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/ultra-long/03_08_22_R941_HG002_rebasecalling-guppy-6.3.7/03_08_22_R941_HG002_4.fq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/ultra-long/03_08_22_R941_HG002_rebasecalling-guppy-6.3.7/03_08_22_R941_HG002_5.fq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/NHGRI_UCSC_panel/HG002/nanopore/ultra-long/03_08_22_R941_HG002_rebasecalling-guppy-6.3.7/03_08_22_R941_HG002_6.fq.gz
```

### Hi-C · ~50× · FASTQ.GZ (paired)

Downsampled Hi-C read pairs from the HPRC panel. Used only in Stage 4 by verkko for phasing. `HIC1` and `HIC2` must be listed in matching order since verkko pairs them positionally.

```bash
# Read 1 files
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S1_R1_001.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S2_R1_001.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S3_R1_001.fastq.gz

# Read 2 files (same order as R1)
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S1_R2_001.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S2_R2_001.fastq.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/hic/downsampled/HG002.HiC_1_S3_R2_001.fastq.gz
```

### Precomputed shortcuts

If you want to enter the pipeline at a later stage, precomputed files are available.

S3 prefix: `s3://human-pangenomics/submissions/09cd8aa1-726c-4cb3-aeea-36e71bab75ff--HG002_hybrid_benchmark/`

| File | Enter at |
|------|----------|
| `r10-hifiasm-correct.ec.ont.fq.gz` | Assembly 4 |
| `assembly.fasta.gz` | Assembly 6 |
| `assembly.haplotype1.fasta.gz`, `assembly.haplotype2.fasta.gz` | Assembly 6 |
| `gqc.tar.gz` | QC · GQC |

---

## Assembly 2 — Define inputs and shell variables

Everything downstream refers to these. Set them once in the shell you will run the rest of the pipeline from, or put them in a small sourced file so a resumed session picks them up.

```bash
HIFI="hg002v1.0.1_hifi_revio_pbmay24.bam"
ONT_R10="hg002v1.0_ont_r10_ul_dorado.bam"
ONT="03_08_22_R941_HG002_4.fq.gz 03_08_22_R941_HG002_5.fq.gz 03_08_22_R941_HG002_6.fq.gz"
HIC1="HG002.HiC_1_S1_R1_001.fastq.gz HG002.HiC_1_S2_R1_001.fastq.gz HG002.HiC_1_S3_R1_001.fastq.gz"
HIC2="HG002.HiC_1_S1_R2_001.fastq.gz HG002.HiC_1_S2_R2_001.fastq.gz HG002.HiC_1_S3_R2_001.fastq.gz"

zipCPUs=8
zipOpt="-l 9 -i"   # -i creates the .gzi index; -I names it
threads=64          # set to the core count on your system
```

Keep `-i` in `zipOpt`. The commands below pass `-I` to name the index file, but `-I` on its own will not create one.

`HIC1` and `HIC2` must list read 1 and read 2 files in the same order, since verkko pairs them positionally.

---

## Assembly 3 — Correction with hifiasm

Only reads over 20 kbp are corrected here, which saves compute time and memory given how deeply covered this dataset is. For a typical dataset, 10 kbp is the recommended threshold. Lower it if your ONT coverage is modest.

### 3.1 — Normalize input formats

```bash
to_fastq() {
  for f in "$@"; do
    case "$f" in
      *.bam)                  samtools fastq "$f" ;;
      *.fastq.gz|*.fq.gz)     zcat "$f" ;;
      *.fastq.bz2|*.fq.bz2)  bzcat "$f" ;;
      *.fastq|*.fq)           cat "$f" ;;
      *)  echo "Warning: unrecognized format for $f" >&2 ;;
    esac
  done
}
```

### 3.2 — Prepare the HiFi input

```bash
to_fastq $HIFI \
  | bgzip -@ $zipCPUs $zipOpt -c -I hifi.input.fastq.gz.gzi - \
  > hifi.input.fastq.gz
```

### 3.3 — Prepare the ONT input

```bash
to_fastq $ONT_R10 \
  | seqtk seq -L 20000 - \
  | bgzip -@ $zipCPUs $zipOpt -c -I r10-hifiasm-correct.input.fastq.gz.gzi - \
  > r10-hifiasm-correct.input.fastq.gz
```

Note this uses `$ONT_R10` only. The optional R9 data is not corrected and enters the pipeline only at the verkko stage.

### 3.4 — Estimate `--hom-cov`

`--hom-cov` is the combined coverage of the ONT and HiFi data going into correction: total bases divided by the genome size (3.1 Gbp for human).

```bash
# With seqkit:
seqkit stats -T hifi.input.fastq.gz r10-hifiasm-correct.input.fastq.gz

# Without seqkit:
for f in hifi.input.fastq.gz r10-hifiasm-correct.input.fastq.gz; do
  zcat "$f" | awk -v f="$f" 'NR%4==2 {n+=length($0)} END {print f, n}'
done
```

Sum the base counts, divide by 3,100,000,000, and round. The value used for this dataset was **185**. This is computed on the length-filtered ONT input from 3.3, not the raw ONT reads.

### 3.5 — Run correction

```bash
hifiasm -e --write-ec \
  --ont r10-hifiasm-correct.input.fastq.gz \
  --hf hifi.input.fastq.gz \
  --hom-cov 185 \
  -o r10-hifiasm-correct.WORKING \
  -t $threads
```

`--ont` supplies the reads to be corrected and `--hf` the reads used to correct them. `-e` with `--write-ec` runs the correction stage and writes corrected reads to disk rather than proceeding to a finished assembly.

The `WORKING` prefix is deliberate — outputs stay under that name until every downstream step succeeds, so an interrupted run cannot leave a truncated file at the final path.

Main output: `r10-hifiasm-correct.WORKING.ec.fq`

### 3.6 — Compress and index

```bash
bgzip -@ $zipCPUs $zipOpt r10-hifiasm-correct.WORKING.ec.fq
samtools faidx r10-hifiasm-correct.WORKING.ec.fq.gz
```

### 3.7 — Split corrected reads by platform

Hifiasm writes HiFi and ONT reads into one corrected file, so they need separating. The split keys on read name prefix: PacBio names begin with the movie ID (`m84039_...`) while ONT names are UUIDs. Since `m` is not a hex character, no UUID can start with it, making `^m` a clean discriminator.

```bash
grep    "^m" r10-hifiasm-correct.WORKING.ec.fq.gz.fai | awk '{print $1}' \
  > r10-hifiasm-correct.hifi.ids
grep -v "^m" r10-hifiasm-correct.WORKING.ec.fq.gz.fai | awk '{print $1}' \
  > r10-hifiasm-correct.ont.ids

wc -l r10-hifiasm-correct.hifi.ids r10-hifiasm-correct.ont.ids
```

!!! warning
    Check the counts before continuing. If either file is empty or the ratio looks wrong, your reads were renamed somewhere upstream and the prefix rule no longer applies. Reads pulled from SRA are the common case, since they arrive with accession-style names.

```bash
seqtk subseq r10-hifiasm-correct.WORKING.ec.fq.gz r10-hifiasm-correct.hifi.ids \
  | bgzip -@ $zipCPUs $zipOpt -c -I r10-hifiasm-correct.ec.hifi.fq.gz.gzi \
  > r10-hifiasm-correct.ec.hifi.fq.gz

seqtk subseq r10-hifiasm-correct.WORKING.ec.fq.gz r10-hifiasm-correct.ont.ids \
  | bgzip -@ $zipCPUs $zipOpt -c -I r10-hifiasm-correct.ec.ont.fq.gz.gzi \
  > r10-hifiasm-correct.ec.ont.fq.gz
```

### 3.8 — Finalize and clean up

```bash
mv r10-hifiasm-correct.WORKING.ec.fq.gz     r10-hifiasm-correct.ec.fq.gz
mv r10-hifiasm-correct.WORKING.ec.fq.gz.gzi r10-hifiasm-correct.ec.fq.gz.gzi
mv r10-hifiasm-correct.WORKING.ec.fq.gz.fai r10-hifiasm-correct.ec.fq.gz.fai
```

Only `r10-hifiasm-correct.ec.ont.fq.gz` is used downstream. Once confirmed intact:

```bash
rm r10-hifiasm-correct.ec.fq.gz*
rm r10-hifiasm-correct.ec.hifi.fq.gz*
rm r10-hifiasm-correct.WORKING.*
```
