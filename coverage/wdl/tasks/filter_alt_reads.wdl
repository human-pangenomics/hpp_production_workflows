version 1.0 

workflow runFilterAltReads {
    call filterAltReads
    output {
        File filteredBam = filterAltReads.filteredBam
    }
}

task filterAltReads {
    input {
        File vcf
        File bam
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_bcftools:latest"
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
        
        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%.bam}
        bcftools view -Ov -m2 -M2 -v snps ~{vcf} > bi_snps.vcf
        mkdir output
        ${FILTER_ALT_READS_BIN} -i ~{bam} -o output/$PREFIX.alt_filtered.bam -v bi_snps.vcf
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File filteredBam = glob("output/*.bam")[0]
    }
}

