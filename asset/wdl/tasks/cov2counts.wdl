version 1.0 

workflow runCov2Counts{
    call cov2counts
}
task cov2counts {
    input {
        File coverageFile
        String sampleName
        String sampleSuffix
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
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

        ${COV2COUNTS_C} -i ~{coverageFile} -o ~{sampleName}.~{sampleSuffix}.counts 
        pigz -p4 -c ~{coverageFile} > ~{sampleName}.~{sampleSuffix}.cov.gz
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File coverageGz = "${sampleName}.${sampleSuffix}.cov.gz"
        File counts = "${sampleName}.${sampleSuffix}.counts"
    }
}

