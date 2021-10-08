version 1.0 

workflow runPhaseReads{
    input {
        File inputBam
        File diploidAssemblyFastaGz
    }
    call sortByName{
         input:
             bamFile = inputBam
    }
    call phaseReads{
         input:
             bamFile = sortByName.outputBam,
             diploidAssemblyFastaGz = diploidAssemblyFastaGz
    }
    call sortByContig{
         input:
             bamFile = phaseReads.outputBam
    }
    output{
        File outputBam = sortByContig.outputBam
        File outputBai = sortByContig.outputBai
        File errLog = phaseReads.errLog
        File outLog = phaseReads.outLog
    }
}

task phaseReads {
    input {
        File bamFile
        File diploidAssemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=500
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
        ${PHASE_READS_BIN} -i ~{bamFile} -f asm.fa -o output/${BAM_PREFIX}.phased.bam 2> err.log > out.log
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

task sortByName {
    input {
        File bamFile
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=500
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
        samtools sort -n -@8 -b ~{bamFile} > output/${BAM_PREFIX}.bam
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
        Int memSize=16
        Int threadCount=8
        Int diskSize=500
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
        samtools sort -@8 -b ~{bamFile} > output/${BAM_PREFIX}.bam
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

