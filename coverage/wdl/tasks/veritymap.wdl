version 1.0 

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t

workflow runVerityMap{
    input {
        File bam
        File bamIndex
        File oneLineHorBed
        File assemblyFastaGz
    }
    call subsetBam{
        input:
            bam = bam,
            bamIndex = bamIndex,
            oneLineBed = oneLineHorBed
    }
    call extractReads_t.extractReads as extractReads {
        input:
            readFile=subsetBam.subsetBam,
            referenceFasta=assemblyFastaGz,
            memSizeGB=4,
            threadCount=4,
            diskSizeGB=2 * ceil(size(subsetBam.subsetBam, "GB")) + 64,
            dockerImage="tpesout/hpp_base:latest"
    }
    call verityMap{
        input:
            assemblyFastaGz = assemblyFastaGz,
            horBed=oneLineHorBed,
            readsFastq = extractReads.extractedRead,
            diskSize = 2 * ceil(size(extractReads.extractedRead, "GB")) + 64,
    }
    output{
        File outputTarGz = verityMap.outputTarGz
    }
}

task subsetBam {
    input {
        File bam
        File bamIndex
        File oneLineBed
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:v0.1"
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
        PREFIX=${FILENAME%%.bam}

        mkdir input
        ln ~{bam} input/${PREFIX}.bam
        ln ~{bamIndex} output/${PREFIX}.bam.bai

        REGION=$(cat ${oneLineBed} | awk '{printf $1":"$2"-"$3}')
        mkdir output
        samtools view -hb input/${PREFIX}.bam ${REGION} > output/${PREFIX}.bam
        samtools index output/${PREFIX}.bam
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File subsetBam = glob("output/*.bam")[0]
        File subsetBamIndex = glob("output/*.bam.bai")[0]
    }
}

task verityMap {
    input {
       	File assemblyFastaGz
        File horBed
        File readsFastq
        String suffix
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="mobinasri/veritymap:latest"
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

        # extract HOR assembly
        FILENAME=$(basename ~{assemblyFastaGz})
        PREFIX=${FILENAME%%.fa*.gz}
        gunzip -c ${FILENAME} > ${PREFIX}.fa
        # make .genome file for 'bedtools getfasta'
        samtools faidx ${PREFIX}.fa
        cat ${PREFIX}.fa.fai | cut -f1-2 > ${PREFIX}.fa.genome

        bedtools getfasta -fi ${PREFIX}.fa -fo ${PREFIX}.~{suffix}.fa -bed ~{horBed} -g ${PREFIX}.fa.genome

        # run VerityMap
        python3 ${VERITY_MAP_PY} --reads ~{readsFastq} -o output -t~{threadCount} -d hifi --careful ${PREFIX}.~{suffix}.fa

        # convert sam to bam, sort and index bam file
        SAM_FILENAME=$(output/*.sam)
        SAM_PREFIX=${SAM_FILENAME%%.sam} 
        samtools view -hb output/${SAM_PREFIX}.sam > output/${SAM_PREFIX}.bam
        samtools sort output/${SAM_PREFIX}.bam > output/${SAM_PREFIX}.sorted.bam
        samtools index output/${SAM_PREFIX}.sorted.bam

        # Rename output folder
        mv output ${PREFIX}.~{suffix}
        tar cvzf ${PREFIX}.~{suffix}.tar.gz ${PREFIX}.~{suffix}
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File outputTarGz = glob("*.tar.gz")[0]
    }
}

