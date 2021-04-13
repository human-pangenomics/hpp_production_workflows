version 1.0

workflow runPstools {
    call pstools
    output {
        File pstoolsOutput = pstools.pstoolsOutput
    }
}

task pstools {
    input {
        Array[File] pairedFastqGzFiles
        File matAssemblyFastaGz
        File patAssemblyFastaGz
        String sampleName
        # runtime configurations
        Int memSize=128
        Int threadCount=16
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_pstools:latest"
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

        mkdir hic_data
        ln ~{sep=" " pairedFastqGzFiles} hic_data/
        pstools phasing_error -t~{threadCount} ~{matAssemblyFastaGz} ~{patAssemblyFastaGz} <(zcat hic_data/*_R1_001.fastq.gz) <(zcat hic_data/*_R2_001.fastq.gz) > ~{sampleName}.pstools.out
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File pstoolsOutput = "${sampleName}.pstools.out"
    }
}


