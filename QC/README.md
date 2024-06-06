# HPRC Assembly QC

*The [Human Pangenome Reference Consortium's](https://humanpangenome.org/) assembly QC pipelines are detailed below.*

This repository includes WDLs run by the HPRC as part of it's production efforts. If you are unable to run WDL workflows each section includes an "equivalent" command which should produce the same or equivalent outputs to the WDL.

## Standard QC 

### Standard QC Outline
![Standard QC](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/docs/imgs/hprc_standard_qc.png?raw=true)

### Tools Used
The following tools are included in the standard_qc pipeline:
* [asmgene](https://github.com/lh3/minimap2)
* [compleasm](https://github.com/huangnengCSU/compleasm)
* [dipcall](https://github.com/lh3/dipcall/tree/v0.2)
* [t2t statistics calculation](https://github.com/biomonika/HPP/blob/main/assembly/wdl/workflows/evaluateHumanAssembly.wdl)
* [yak](https://github.com/lh3/yak)
* [merqury](https://github.com/marbl/merqury)
* [paftools' misjoin check](https://github.com/lh3/minimap2/blob/67dd906a80988dddacc8c551623fdc75b0c12dd2/misc/paftools.js#L2605-L2719)

## Self Alignment Based QC
![Alignment Based QC](https://github.com/human-pangenomics/hpp_production_workflows/blob/master/docs/imgs/hprc_alignment_based_qc.png?raw=true)

### Tools Used

The following tools are included in the alignment_based_qc pipeline:
* [flagger](https://github.com/mobinasri/flagger/tree/main)
* [nucfreq](https://github.com/mrvollger/NucFreq/tree/master)

