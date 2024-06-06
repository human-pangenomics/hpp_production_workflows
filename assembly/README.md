# HPRC Assembly Production

*The [Human Pangenome Reference Consortium's](https://humanpangenome.org/) assembly production pipeline is detailed below.*

![Assembly Process](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/docs/imgs/hprc_assembly_steps.jpeg?raw=true)

This repository includes WDLs run by the HPRC as part of it's production efforts. If you are unable to run WDL workflows each section includes an "equivalent" command which should produce the same or equivalent outputs to the WDL.

## Assembly Workflow

Two workflows are used to create phased assemblies:
* trio_hifiasm_assembly_cutadapt_multistep.wdl
* hic_hifiasm_assembly_cutadapt_multistep.wdl

If you would like to see example input jsons for these WDLs they can be found in the `example_inputs` folder:
* wdl/example_inputs/trio_hifiasm_assembly_cutadapt_multistep_example_inputs.json
* wdl/example_inputs/hic_hifiasm_assembly_cutadapt_multistep_example_inputs.json

**The workflows have three main steps**
* Cutadapt: trim/filter adapters from HiFi reads
* Yak (Trio only): build kmer database from parental Illumina reads
* Hifiasm: create assemblies using either trio or Hi-C phasing
* GroupXY (Hi-C only): correct sex chromosome assignment in male samples


### Cutadapt

PacBio Hifi adapters are rare in reads, but can be incorporated into the assembly if they are present. To remove them (and avoid having to hard mask out adapter sequences from the assembly) [Cutadapt](https://github.com/marcelm/cutadapt/) is used to discard any reads containing adapter sequence.

<details>
<summary>Equivalent Command</summary>
<br>

```
cutadapt \
    -b "AAAAAAAAAAAAAAAAAATTAACGGAGGAGGAGGA;min_overlap=35" \
    -b "ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT;min_overlap=45" \
    --discard-trimmed \
    -o ${PREFIX}.filt.fastq.gz \
    ${readFastqGz} \
    -j ${threadCount} \
    --revcomp \
    -e 0.05
```

</details>


### yak count

When parental Illumina data is available [yak](https://github.com/lh3/yak) is used to create kmer databases that are passed to hifiasm for trio binning the assemblies.

<details>
<summary>Equivalent Command</summary>
<br>

```
bloomSize=37
readFiles=$(IFS=' '; echo "${readFiles[*]}") 

yak count \
    -t${threadCount} \
    -b${bloomSize} \
    -o ${sampleName}.yak \
    <(cat ${readFiles}) <(cat ${readFiles})
```

</details>


### Hifiasm

Assemblies are produced from PacBio HiFi and ONT Ultralong reads using [Hifiasm](https://github.com/chhylp123/hifiasm). Either parental Illumina data or Hi-C data is used for phasing to produce phased, T2T (or near T2T) contigs and scaffolds.

A few parameters are added to the hifiasm call to improve performance:
* `--telo-m CCCTAA`: helps to produce more telomeric sequence at the ends of contigs/scaffolds
* `--dual-scaf`: scaffold together contigs based on the graph's structure (graphold)
* `--ul-cut`: filter out ultralong reads below cutoff

<details>
<summary>Equivalent Trio Command</summary>
<br>

```
minOntReadLength=50000
childReadsUL=$(IFS=,; echo "${childReadsUL[*]}")
childReadsHiC1=$(IFS=' '; echo "${childReadsHiC1[*]}") 
childReadsHiC2=$(IFS=' '; echo "${childReadsHiC2[*]}") 
childReadsHiFi=$(IFS=' '; echo "${childReadsHiFi[*]}") 

hifiasm \
    --telo-m CCCTAA \
    --dual-scaf \
    --ul-cut ${minOntReadLength} \
    -t${threadCount} \    
    -o ${childID} \
    --ul "${childReadsUL}" \
    --hom-cov ${homCov} \
    -1 "${paternalYak}" \
    -2 "${maternalYak}" \
    "${childReadsHiFi}"
```

</details>

<details>
<summary>Equivalent HiC Command</summary>
<br>

```
minOntReadLength=50000
childReadsUL=$(IFS=,; echo "${childReadsUL[*]}")
childReadsHiC1=$(IFS=' '; echo "${childReadsHiC1[*]}") 
childReadsHiC2=$(IFS=' '; echo "${childReadsHiC2[*]}") 
childReadsHiFi=$(IFS=' '; echo "${childReadsHiFi[*]}") 

hifiasm \
    --telo-m CCCTAA \
    --dual-scaf \
    --ul-cut ${minOntReadLength} \
    -t${threadCount} \    
    -o ${childID} \
    --ul "${childReadsUL}" \
    --hom-cov ${homCov} \
    --h1 "${childReadsHiC1}" \
    --h2 "${childReadsHiC2}"  \
    "${childReadsHiFi}"
```

</details>

### GroupXY (Hi-C Only)

Hi-C phased assemblies produced by hifiasm place contigs from chromosome Y in both haplotypes in male samples. To correct sex chromosome assignment Yak's [groupxy](https://github.com/lh3/yak) command is used on male assembles.


<details>
<summary>Equivalent Command</summary>
<br>

```
yak sexchr \
    -K2g \
    -t16 \
    ${chrY_no_par_yak} \
    ${chrX_no_par_yak} \
    ${par_yak} \
    ${hap1_gz} \
    ${hap2_gz} \
    > cnt.txt

groupxy.pl \
    cnt.txt \
    | awk '$4==1' | cut -f2 \
        | seqtk subseq -l60 \
        <(zcat ${hap1_gz} ${hap2_gz}) - \
        | pigz \
        > ${childID}.hap1.corrected.fa.gz

groupxy.pl \
    cnt.txt | \
    awk '$4==2' | cut -f2 \
        | seqtk subseq -l60 \
        <(zcat ${hap1_gz} ${hap2_gz}) - \
        | pigz \
        > ${childID}.hap2.corrected.fa.gz
```

</details>

