version 1.0 

workflow runBed2Fasta{
    input {
        String sampleName
        String sampleSuffix
        File assemblyFastaGz
        File bed
    }
    call extendBedAndExtractSequence as bedTask{
        input:
            sampleName = "~{sampleName}.~{sampleSuffix}",
            assemblyFastaGz = assemblyFastaGz,
            bed = bed
    }
    output {
        File extendedBlocksBED = bedTask.extendedBlocksBED
        File extractedBlocksFastaGz = bedTask.extractedBlocksFastaGz
    }
}
task extendBedAndExtractSequence {
    input {
        String sampleName
        File assemblyFastaGz
        File bed
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

        gunzip -c ~{assemblyFastaGz} > ~{sampleName}.fa

        BED_FILENAME=$(basename ~{bed})
        ln -s ~{bed} ${BED_FILENAME}
        $BLOCKS_EXTENSION_V1_BASH ${BED_FILENAME} ~{sampleName}

        bedtools getfasta -fi ~{sampleName}.fa -bed *.extended.bed > ~{sampleName}.extended.fa
        gzip ~{sampleName}.extended.fa
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File extendedBlocksBED = glob("*.extended.bed")[0]
        File extractedBlocksFastaGz = "~{sampleName}.extended.fa.gz"
    }
}

