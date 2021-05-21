version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t

workflow RunDownSampling{
    input {
        File readFile
        Int downSampledCoverage
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
    }

    call extractReads_t.extractReads as extractReads {
        input:
            readFile=readFile,
            referenceFasta=referenceFasta,
            memSizeGB=4,
            threadCount=4,
            diskSizeGB=fileExtractionDiskSizeGB,
            dockerImage="tpesout/hpp_base:latest"
    }
    call seqtkTask {
        input:
            readFile=extractReads.extractedRead,
            downSampledCoverage=downSampledCoverage,
            diskSize=2 * extractReads.fileSizeGB
    }
    output {
        File downSampledFastq=seqtkTask.downSampledFastq
        File totalBasesStats=seqtkTask.totalBasesStats
    }
}

task seqtkTask {
    input{
        File readFile
        Int downSampledCoverage
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=256
        Int preemptible=2
        String dockerImage="quay.io/masri2019/hpp_seqtk:latest"
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

        mkdir downsampled
        FILE_NAME=$(basename ~{readFile})
        PREFIX=${FILE_NAME%.*}
    
        echo "Total number of bases before down-sampling:" > totalBases.txt
        cat ~{readFile} | paste - - - - | cut -f 2 | tr -d '\n' | wc -c >> totalBases.txt
        SAMPLING_RATE=$(sed -n 2p totalBases.txt | awk -v downsampledCoverage="~{downSampledCoverage}" '{coverage=$1/3.1e9;  printf "%.2f\n",downsampledCoverage / coverage}') 
        seqtk sample ~{readFile} ${SAMPLING_RATE} > downsampled/${PREFIX}.~{downSampledCoverage}X.fq
        echo "Total number of bases after down-sampling:" >> totalBases.txt
        cat downsampled/${PREFIX}.~{downSampledCoverage}X.fq | paste - - - - | cut -f 2 | tr -d '\n' | wc -c >> totalBases.txt
        
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
        File totalBasesStats = "totalBases.txt"
        File downSampledFastq = glob("downsampled/*")[0]
    }
}
