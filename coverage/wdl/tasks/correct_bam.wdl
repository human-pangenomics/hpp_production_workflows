version 1.0 

workflow runCorrectBam{
    call correctBam
    output {
        File correctedBam = correctBam.correctedBam
    }
}

task correctBam {
    input {
        File? phasingLogText
        File? mapqTableText
        File bam
        String? options
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=512
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
        
        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%.bam}

        mkdir output
        OPTIONS="~{options}"

        if [ -n "~{phasingLogText}" ]
        then
            OPTIONS="${OPTIONS} --phasingLog ~{phasingLogText}"
        fi

        if [ -n "~{mapqTableText}" ]
        then
            OPTIONS="${OPTIONS} --mapqTable ~{mapqTableText}"
        fi
        
        ${CORRECT_BAM_BIN} ${OPTIONS} -i ~{bam} -o output/$PREFIX.corrected.bam
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File correctedBam = glob("output/*.bam")[0]
    }
}

