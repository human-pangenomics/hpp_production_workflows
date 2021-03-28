# HPP Production Workflows

This repository holds WDL workflows and Docker build scripts for 
production workflows used by the [Human Pangenome Project](https://humanpangenome.org/).

### Organization
The repo is split into assembly and QC workflows with the following organization:

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

------------------


## Assembly

Assemblies are produced with Hifiasm from reads filtered with HiFiAdapterFilt:
* Filter HiFi Data w/ [HiFiAdapterFilt](https://github.com/sheinasim/HiFiAdapterFilt/tree/master/DB)
* Run [Hifiasm](https://github.com/chhylp123/hifiasm)

------------------


## QC

### Automated Assembly QC

Assembly QC is run with one of two workflows:
* standard_qc.wdl
* standard_qc_no_qv.wdl (for samples with no child Illumina data)

The following tools are included:
* [asmgene](https://github.com/lh3/minimap2)
* [dipcall v0.1](https://github.com/lh3/dipcall/tree/v0.1)
* [dipcall v0.2](https://github.com/lh3/dipcall/tree/v0.2)
* [merqury](https://github.com/marbl/merqury)
* [QUAST](https://sourceforge.net/projects/quast/files/)
* [YAK](https://github.com/lh3/yak)

### Miscelaneous QC

In addition to the automated QC of assemblies, the following tasks are included:
* read_stats.wdl: outputs read length statistics for bam/fastq files
* mask_assembly.wdl: masks assemblies in regions found in bed file
* dropFastaContigs.wdl: drops contigs from assemblies based on text file list of contigs to drop
------------------ 


