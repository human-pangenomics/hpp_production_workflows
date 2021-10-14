version 1.0

workflow runGetPhasedReadIds{
    call getPhasedReadIds
    output{
        File hap1IdsText = getPhasedReadIds.hap1IdsText
        File hap2IdsText = getPhasedReadIds.hap2IdsText
        File ambiguousIdsText = getPhasedReadIds.ambiguousIdsText
    }  
}
task getPhasedReadIds{
    input{
        File bamFile
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=256
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

        FILENAME=$(basename ~{bamFile})
        PREFIX=${FILENAME%.bam}

        mkdir output
        samtools view -F4 -F256 ~{bamFile} | awk '{print $1"\t"$3}' > read_ids_contig.txt 
        cat read_ids_contig.txt | awk '$2 ~ /.*#1#.*/ {print $1}' | sort -u > hap1.txt
        cat read_ids_contig.txt | awk '$2 ~ /.*#2#.*/ {print $1}' | sort -u > hap2.txt
        comm -23 hap1.txt hap2.txt > output/${PREFIX}.hap1_read_ids.txt
        comm -13 hap1.txt hap2.txt > output/${PREFIX}.hap2_read_ids.txt
        comm -12 hap1.txt hap2.txt > output/${PREFIX}.ambiguous_read_ids.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File hap1IdsText = glob("output/*.hap1_read_ids.txt")[0]
        File hap2IdsText = glob("output/*.hap2_read_ids.txt")[0]
        File ambiguousIdsText = glob("output/*.ambiguous_read_ids.txt")[0]
    }
}

