version 1.0 

workflow runCalmd{
    call calmd
    output{
        File mdBamFile = calmd.mdBamFile
    }
}
task calmd {
    input {
        File bamFile
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=500
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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
        
        FILENAME=$(basename ~{bamFile})
        PREFIX=${FILENAME%.bam}
        samtools calmd -@8 -b ~{bamFile} > ${PREFIX}.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mdBamFile = glob("*.bam")[0]
    }
}

