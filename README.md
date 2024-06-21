# HPP Production Workflows

This repository holds WDL workflows and Docker build scripts for
production workflows for data QC, assembly generation, and assembly QC used by the [Human Pangenome Reference Consortium](https://humanpangenome.org/).

All WDLs and containers created in this repository are licensed under the MIT license. The underlying tools (that the WDLs and containers run) are likely covered under one or more Free and Open Source Software licenses, but we cannot make any guarantees to that fact.

------------------

## Repository Organization
Workflows are split across data_processing, assembly, and (assembly) QC folders; each with the following folder structure:

```
 ── docker/
    └── toolName/
        └── Dockerfile
        └── Makefile
        └── scripts/
            └── toolName/
                └── scriptName.py
 ── wdl/
    └── tasks/
    │   └── taskName.wdl
    └── workflows/
        └── workFlowName.wdl
```

The root level of the data_processing, assembly, and (assembly) QC folders each contain a readme that provides details about the workflows and how to use them. Summaries of the workflows in each area are below.

------------------

# Workflow Types
## Data Processing

The HPRC produces HiFi, ONT, and Illumina Hi-C data. Each data type has a workflow to check data files to ensure they pass QC.
* HiFi QC Workflow   
   * Check for file-level sample swaps with [NTSM](https://github.com/JustinChu/ntsm)
   * Calculate coverage (Gbp) and insert(N50) metrics from fastqs/bams using [in-house tooling](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/QC/docker/read_stats/scripts/fai_read_stats/fai_read_stats.py)
   * Check for methylation and kinetics tags (in progress)
* ONT QC Workflow   
   * Check for file-level sample swaps with [NTSM](https://github.com/JustinChu/ntsm)
   * Calculate coverage (Gbp) and insert(N50) metrics from summary files using [in-house tooling](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/data_processing/docker/ont_summary_stats/scripts/calculate_summary_stats.py)
* Hi-C QC Workflow   
   * Check for file-level sample swaps with [NTSM](https://github.com/JustinChu/ntsm)
   * Calculate total bases for the data file

## Assembly

Assemblies are produced with one of two [Hifiasm](https://github.com/chhylp123/hifiasm) workflows using HiFi and ONT ultralong reads with phasing by either Illumina Hi-C or parental Illumina data for the Hi-C and trio workflows, respectively. The major steps included in the assembly workflows are:
* [Yak](https://github.com/lh3/yak) for creation of kmer databases for trio phased assemblies
* [Cutadapt](https://github.com/marcelm/cutadapt) for adapter filtering of HiFi reads
* Run Hifiasm with HiFi and ONT ultralong and trio or Hi-C phasing
* [Yak](https://github.com/lh3/yak) for sex chromosome assignment in Hi-C phased assemblies

In addition to the Hifiasm workflows there is an assembly cleanup workflow which:
* Removes contamination with [NCBI's FCS](https://github.com/ncbi/fcs)
* Removes mitochondrial contigs
* Runs [MitoHiFi](https://github.com/marcelauliano/MitoHiFi) to assemble mitochondrial contigs
* Assigns chromosome labels to fasta headers of T2T contigs/scaffolds

## Polishing

Assemblies are polished using a custom pipeline based around [DeepPolisher](https://github.com/google/deeppolisher). The major steps in the HPRC assembly polishing pipeline are:
* Alignment of all HiFi reads to the diploid assembly using [minimap2](https://github.com/lh3/minimap2)
* Alignment of all ONT UL reads > 100kb separately to each haplotype assembly using [minimap2](https://github.com/lh3/minimap2)
* [PHARAOH pipeline](https://github.com/miramastoras/PHARAOH). PHARAOH ensures optimal HiFi read phasing, by leveraging ONT UL information to assign HiFi reads to the correct haplotype in stretches of homozygosity longer than 20kb.
* [DeepPolisher](https://github.com/google/deeppolisher) is an encoder-only transformer model which is run on the PHARAOH-corrected HiFi alignments, to predict polishing edits in the assemblies.

## QC

### Automated Assembly QC

Assembly QC is broken down into two types:
* standard_qc: these tools are relatively fast to run and provide insight into the completeness, correctness, and contiguity of the assemblies.
* alignment_based_qc: these tools rely on long read alignment of a sample's reads to it's own assembly. The alignments are then used to identify unexpected variation that indicates missassembly.

The following tools are included in the standard_qc pipeline:
* [asmgene](https://github.com/lh3/minimap2)
* [compleasm](https://github.com/huangnengCSU/compleasm)
* [dipcall](https://github.com/lh3/dipcall/tree/v0.2)
* [t2t statistics calculation](https://github.com/biomonika/HPP/blob/main/assembly/wdl/workflows/evaluateHumanAssembly.wdl)
* [yak](https://github.com/lh3/yak)
* [merqury](https://github.com/marbl/merqury)
* [paftools' misjoin check](https://github.com/lh3/minimap2/blob/67dd906a80988dddacc8c551623fdc75b0c12dd2/misc/paftools.js#L2605-L2719)

The following tools are included in the alignment_based_qc pipeline:
* [flagger](https://github.com/mobinasri/flagger/tree/main)
* [nucfreq](https://github.com/mrvollger/NucFreq/tree/master)

------------------

# Running WDLs
If you haven't run a WDL before, there are good resources online to get started. You first need to choose a way to run WDLs. Below are a few options:
* [Terra](https://app.terra.bio/): An online platform with a GUI which can run workflows in Google Cloud or Microsoft Azure.
* [Cromwell](https://cromwell.readthedocs.io/en/latest/tutorials/FiveMinuteIntro/): A workflow runner from the BROAD with support for Slurm, local compute, and multiple cloud platforms.
* [Toil](https://toil.readthedocs.io/en/master/wdl/running.html): A workflow runner from UCSC with support for Slurm, local compute, and multiple cloud platforms.

## Running with Cromwell  

Before starting, read the [Cromwell 5 minute intro](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/).

Once you've done that, download the latest version of cromwell and make it executable. (Replace XY with newest version number)  
```
wget https://github.com/broadinstitute/cromwell/releases/download/86/cromwell-XY.jar
chmod +x cromwell-XY.jar
```

And run your WDL:  
```
java -jar cromwell-XY.jar run \
   /path/to/my_workflow.wdl \
   -i my_workflow_inputs.json \
   > run_log.txt
```

### Input files

Each workflow requires an input json. You can create a template using womtool:

```
java -jar womtool-XY.jar \
    inputs \
    /path/to/my_workflow.wdl \
    > my_workflow_inputs.json
```
------------------
