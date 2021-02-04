version 1.0

workflow runFilterHiFiAdapter{
    call filterHiFiAdapter
}


task filterHiFiAdapter {
    input{
        File readFastq
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=16
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/masri2019/hpp_hifi_adapter_filt:latest"
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

        mkdir data
        cd data
        FILENAME=$(basename -- "~{readFastq}")
        PREFIX="${FILENAME%.*}"
        ln ~{readFastq} ${PREFIX}.fastq
        bash ${HIFI_ADAPTER_FILTER_BASH} -t ~{threadCount}
        OUTPUTSIZE=`du -s -BG *.filt.fastq | sed 's/G.*//'`
        echo $OUTPUTSIZE > outputsize
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File filteredReadFastq = glob("data/*.filt.fastq")[0]
        File blastout = glob("data/*.contaminant.blastout")[0]
        File blocklist = glob("data/*.blocklist")[0] 
        Int fileSizeGB = read_int("data/outputsize") 
    }
}

