version 1.0 

import "bam2paf.wdl" as bam2paf_t

workflow runProjectBlocks {
    input {
        String assemblyFastaGz
        File asm2refBam
        File refBlocksBed
        String suffix
    }
    call bam2paf_t.bam2paf {
       input:
           bamFile = asm2refBam,
           minMAPQ = 0,
           primaryOnly = "yes"
    }
    call projectRef2Asm {
        input:
            refBlocksBed = refBlocksBed,
            asm2refPaf = bam2paf.pafFile,
            sampleName = basename("${assemblyFastaGz}", ".fa.gz"),
            suffix = suffix
    }
    output {
       File projectionBed = projectRef2Asm.projectionBed
    }
}

task projectRef2Asm {
    input {
        File refBlocksBed
        File asm2refPaf
        String sampleName
        String suffix
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        python3 ${PROJECT_BLOCKS_PY} --mode 'ref2asm' --paf ~{asm2refPaf} --blocks ~{refBlocksBed} --outputProjectable projectable.bed --outputProjection ~{sampleName}.~{suffix}.bed
        mkdir output
        bedtools sort -i ~{sampleName}.~{suffix}.bed | bedtools merge -i - > output/~{sampleName}.~{suffix}.bed

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File projectionBed = glob("output/*.bed")[0]
    }
}
