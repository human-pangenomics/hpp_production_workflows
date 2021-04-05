version 1.0 

workflow runMergeAssetOutputsLongRead{
    input {
        String sampleName
        String sampleSuffix
        File assemblyFastaGz
        File assetOutputsTarGz
        File cenSatRegionsBed
    }
    call mergeAssetOutputs as merge{
        input:
            sampleName = "~{sampleName}.~{sampleSuffix}",
            assemblyFastaGz = assemblyFastaGz,
            assetOutputsTarGz = assetOutputsTarGz,
            cenSatRegionsBed = cenSatRegionsBed
    }
    output {
        File blocksStatsTEXT = merge.blocksStatsTEXT
        File assetLowSupportBED = merge.assetLowSupportBED
    }
}

task mergeAssetOutputs {
    input {
        String sampleName
        File assemblyFastaGz
        File assetOutputsTarGz
        File cenSatRegionsBed
        # runtime configurations
        Int memSize=4
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

        tar -xzf ~{assetOutputsTarGz} --strip-components 1

        gunzip -c ~{assemblyFastaGz} > ~{sampleName}.fa

        $MERGE_BLOCKS_BASH_2 ~{sampleName}
        $BLOCKS_STAT_V2_BASH ~{sampleName} ~{cenSatRegionsBed} > ~{sampleName}.blocks.stats.txt

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File blocksStatsTEXT = "~{sampleName}.blocks.stats.txt"
        File assetLowSupportBED = "~{sampleName}.low_support.trim1k.bed"
    }
}

