version 1.0

import "bam2paf.wdl" as bam2paf_t
import "project_blocks.wdl" as projectBlocks_t

workflow runProjectBlocksMale {
    input {
        String assemblyFastaGz
        File asm2refBam
        File nonCntr_nonSex_nonMito_bed
        File cntr_nonSex_bed
        File chrX_PAR_bed
        File chrX_cntr_bed
        File chrX_nonCntr_nonPAR_bed
        File chrY_PAR_bed
        File chrY_nonPAR_bed
    }
    call bam2paf_t.bam2paf {
       input:
           bamFile = asm2refBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call projectBlocks_t.projectRef2Asm as project1 {
        input:
            refBlocksBed = nonCntr_nonSex_nonMito_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "nonCntr_nonSex_nonMito"
    }
    call projectBlocks_t.projectRef2Asm as project2 {
        input:
            refBlocksBed = cntr_nonSex_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "cntr_nonSex"
    }
    call projectBlocks_t.projectRef2Asm as project3 {
        input:
            refBlocksBed = chrX_PAR_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "chrX_PAR_bed"
    }
    call projectBlocks_t.projectRef2Asm as project4 {
        input:
            refBlocksBed = chrX_cntr_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "chrX_cntr"
    }
    call projectBlocks_t.projectRef2Asm as project5 {
        input:
            refBlocksBed = chrX_nonCntr_nonPAR_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "chrX_nonCntr_nonPAR"
    }
    call projectBlocks_t.projectRef2Asm as project6 {
        input:
            refBlocksBed = chrY_PAR_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "chrY_PAR"
    }
    call projectBlocks_t.projectRef2Asm as project7 {
        input:
            refBlocksBed = chrY_nonPAR_bed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = "chrY_nonPAR"
    }
    output {
       File projection_nonCntr_nonSex_nonMito_bed = project1.projectionBed
       File projection_cntr_nonSex_bed = project2.projectionBed
       File projection_chrX_PAR_bed = project3.projectionBed
       File projection_chrX_cntr_bed = project4.projectionBed
       File projection_chrX_nonCntr_nonPAR_bed = project5.projectionBed
       File projection_chrY_PAR_bed = project6.projectionBed
       File projection_chrY_nonPAR_bed = project7.projectionBed
    }
}
