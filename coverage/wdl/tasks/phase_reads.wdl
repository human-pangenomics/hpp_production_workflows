version 1.0 

import "merge_bams.wdl" as mergeBams_t

workflow runPhaseReads{
    input {
        File inputBam
        File diploidAssemblyFastaGz
        String sampleName
        String sampleSuffix
    }
    call sortByName{
         input:
             bamFile = inputBam
    }
    call splitByName{
         input:
             bamFile = sortByName.outputBam
    }
    scatter (splitBam in splitByName.splitBams) { 
        call phaseReads{
            input:
                bamFile = splitBam,
                diploidAssemblyFastaGz = diploidAssemblyFastaGz
        }
        call sortByContig{
            input:
                bamFile = phaseReads.outputBam
        }
    }
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            sortedBamFiles = sortByContig.outputBam
    }
    output{
        File outputBam = mergeBams.mergedBam
        File outputBai = mergeBams.mergedBai
        Array[File] errLogs = phaseReads.errLog
        Array[File] outLogs = phaseReads.outLog
    }
}

task phaseReads {
    input {
        File bamFile
        File diploidAssemblyFastaGz
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        gunzip -c ~{diploidAssemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        mkdir output
        ${PHASE_READS_BIN} -q 1 -d 1e-5 -e 0.1 -b 10 -i ~{bamFile} -f asm.fa -o output/${BAM_PREFIX}.phased.bam 2> err.log > out.log
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.phased.bam")[0]
        File errLog = "err.log"
        File outLog = "out.log"
    }
}

task splitByName {
    input {
        File bamFile
        Int NReadsPerBam = 500000
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        ${SPLIT_BAM_BY_READNAME_BIN} -i ~{bamFile} -o output -n ~{NReadsPerBam}
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        Array[File] splitBams = glob("output/*.bam")
    }
}


task sortByName {
    input {
        File bamFile
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_base:latest"
        Int preemptible=2
        String zones="us-west2-a"
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
        
        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}
       
        mkdir output
        samtools sort -n -@8 ~{bamFile} > output/${BAM_PREFIX}.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.bam")[0]
    }
}

task sortByContig {
    input {
        File bamFile
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_base:latest"
        Int preemptible=2
        String zones="us-west2-a"
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

        BAM_FILENAME=$(basename ~{bamFile})
        BAM_PREFIX=${BAM_FILENAME%.bam}

        mkdir output
        samtools sort -@~{threadCount} ~{bamFile} > output/${BAM_PREFIX}.bam
        samtools index output/${BAM_PREFIX}.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones : zones
    }
    output {
        File outputBam = glob("output/*.bam")[0]
        File outputBai = glob("output/*.bai")[0]
    }
}

