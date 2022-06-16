version 1.0 

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "correct_bam.wdl" as correct_bam_t

workflow runVerityMap{
    input {
        File bam
        File bamIndex
        Array[File] oneLineHorBedArray
        Array[File] suffixArray
        File assemblyFastaGz
        File? phasingLogText
        Int minReadLength = 5000
        Int minAlignmentLength = 2000
        Float maxDivergence = 0.01
    }

    ## Correct the bam file by swapping pri/sec tags for the wrongly phased reads
    call correct_bam_t.correctBam {
        input:
            bam = bam,
            phasingLogText = phasingLogText,
            suffix = "corrected",
            options = "--primaryOnly --minReadLen ${minReadLength} --minAlignment ${minAlignmentLength} --maxDiv ${maxDivergence}",
            flagRemoveSupplementary = true,
            flagRemoveMultiplePrimary = true,
            diskSize = ceil(size(bam, "GB")) * 2 + 64
    }

    scatter (bedAndSuffix in zip(oneLineHorBedArray, suffixArray)) {
        File bed = bedAndSuffix.left
        File suffix = bedAndSuffix.right

        # subset the bam file to include only the alignments to the HOR array
        call subsetBam{
            input:
                bam = correctBam.correctedBam,
                bamIndex = correctBam.correctedBamIndex,
                oneLineBed = bed
        }

        # Extract reads in fastq format
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=subsetBam.subsetBam,
                referenceFasta=assemblyFastaGz,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=2 * ceil(size(subsetBam.subsetBam, "GB")) + 64,
                dockerImage="tpesout/hpp_base:latest"
        }

        # Run VerityMap to evaluate the HOR array
        call verityMap{
            input:
                assemblyFastaGz = assemblyFastaGz,
                horBed = bed,
                readsFastq = extractReads.extractedRead,
                suffix = suffix,
                diskSize = 2 * ceil(size(extractReads.extractedRead, "GB")) + 64,
        }
    }
    output{
        Array[File] outputTarGzArray = verityMap.outputTarGz
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
        String dockerImage="mobinasri/bio_base:latest"
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
        ln ~{bamIndex} input/${PREFIX}.bam.bai

        REGION=$(cat ~{oneLineBed} | awk '{printf $1":"$2"-"$3}')
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
        gunzip -c ~{assemblyFastaGz} > ${PREFIX}.fa
        samtools faidx ${PREFIX}.fa

        bedtools getfasta -fi ${PREFIX}.fa -fo ${PREFIX}.~{suffix}.fa -bed ~{horBed}

        # run VerityMap
        python3 ${VERITY_MAP_PY} --reads ~{readsFastq} -o output -t~{threadCount} -d hifi --careful ${PREFIX}.~{suffix}.fa

        # convert sam to bam, sort and index bam file
        SAM_PATH=$(ls output/*.sam)
        SAM_FILENAME=$(basename ${SAM_PATH})
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

