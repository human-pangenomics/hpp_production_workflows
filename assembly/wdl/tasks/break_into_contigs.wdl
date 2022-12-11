version 1.0 

workflow runBreakIntoContigs{
    call breakIntoContigs
}
task breakIntoContigs {
    input {
        File assemblyFaGz
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_hifiasm:0.18.2"
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
        
        FILENAME=$(basename ~{assemblyFaGz})
        PREFIX=${FILENAME%%.fa.gz}

        mkdir output
        python3 /home/scripts/break_into_contigs.py --inputFasta <(zcat ~{assemblyFaGz}) | pigz -p8 - > output/${PREFIX}.contigs.fa.gz
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File assemblyContigsFaGz = glob("output/*.gz")[0]
    }
}

