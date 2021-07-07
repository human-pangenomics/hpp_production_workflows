version 1.0

import "bam2paf.wdl" as bam2paf_t
import "project_blocks.wdl" as projectBlocks_t

workflow runProjectBlocksFemale {
    input {
        String assemblyFastaGz
        File asm2refBam
        File nonCntr_nonMito_bed
        File cntr_bed
    }
    call bam2paf_t.bam2paf {
       input:
           bamFile = asm2refBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call projectBlocks_t.projectRef2Asm as project1 {
        input:
            refBlocksBed = nonCntr_nonMito_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "nonCntr_nonMito"
    }
    call projectBlocks_t.projectRef2Asm as project2 {
        input:
            refBlocksBed = cntr_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "cntr"
    }
    output {
       File projection_nonCntr_nonMito_bed = project1.projectionBed
       File projection_cntr_bed = project2.projectionBed
    }
}
