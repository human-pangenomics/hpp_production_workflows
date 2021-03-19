version 1.0 

workflow runMapBlocks{
    call mapBlocks
    output {
        File mappedBlocksBed = mapBlocks.mappedBlocksBed
    }
}
task mapBlocks {
    input {
        File blocksBed
        File alignmentBam
        String suffix
        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=32
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
        
        FILENAME=$(basename ~{blocksBed})
        PREFIX=${FILENAME%.bed}
        samtools view -F256 -F4 -q20 ~{alignmentBam} | cut -f1-6 > no_seq.sam
        mkdir output
        python3 $MAP_BLOCKS_PY --sam no_seq.sam --bed ~{blocksBed} --output output/${PREFIX}.~{suffix}.bed
        bedtools sort -i output/${PREFIX}.~{suffix}.bed | bedtools merge -i - > output/${PREFIX}.~{suffix}.merged.bed
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mappedBlocksBed = glob("output/*.merged.bed")[0]
    }
}

