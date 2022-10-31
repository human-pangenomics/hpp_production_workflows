version 1.0

import "../../../coverage/wdl/tasks/merge_bams.wdl" as merge_bams_t
import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t

workflow bwaPairedAlignment{
    input {
        String sampleName
        Array[File] readFiles_1
        Array[File] readFiles_2
        File assembly
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 256
        Int preemptible=2
        String zones
    }

    ## build bwa index files for the assembly
    call buildBwaIndex{
        input:
            assembly = assembly,
            # runtime configurations
            preemptible = preemptible,
            zones = zones
     }

    scatter (pairedReadFile in zip(readFiles_1, readFiles_2)){
        call extractReads_t.extractReads as extractReads_1 {
            input:
                readFile=pairedReadFile.left,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage="tpesout/hpp_base:latest"
        }
        call extractReads_t.extractReads as extractReads_2 {
            input:
                readFile=pairedReadFile.right,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage="tpesout/hpp_base:latest"
        }
        call pairedAlignment{
            input:
                referenceFasta = referenceFasta,
                pairedFastq_1 = extractReads_1.extractedRead,
                pairedFastq_2 = extractReads_2.extractedRead,
                indexTar = buildBwaIndex.indexTar,
                diskSize = extractReads_1.fileSizeGB * 6 + 32,
                preemptible = preemptible,
                zones = zones
        }
    }

    call arithmetic_t.sum as bamSize {
        input:
            integers=pairedAlignment.fileSizeGB
    }

    ## merge the bam files
    call merge_bams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sortedBamFiles = pairedAlignment.sortedBamFile,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5),
            preemptible = preemptible,
            zones = zones
    }

    output{
        File sortedBamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }
 
}

task buildBwaIndex{
    input{
        File assembly
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_bwa:latest"
        Int preemptible=2
        String zones
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
        
        if [[ ~{assembly} =~ .*f(ast)?a\.gz$ ]] ; then    
            zcat ~{assembly} > asm.fa
        elif [[ ~{assembly} =~ .*f(ast)?a$ ]] ; then
            ln ~{assembly} asm.fa
        else
             echo "UNSUPPORTED READ FORMAT (expect .fa .fasta .fa.gz .fasta.gz): $(basename ~{assembly})"
             exit 1
        fi

        # build bwa index for the given assembly
        bwa index asm.fa
        mkdir index
        mv asm.* index/
        tar -cf index.tar index
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
        File indexTar = "index.tar"
    }
}

task pairedAlignment{
    input{
        File pairedFastq_1
        File pairedFastq_2
        File indexTar
        String pairedSuffix_1
        File? referenceFasta
        String bwaParams = "-SP -B10"
        # runtime configurations
        Int memSize=64
        Int threadCount=32
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_bwa:latest"
        Int preemptible=2
        String zones
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

        # extract the previously generated bwa index
        tar -xf ~{indexTar} --strip-components 1

        fileBasename=$(basename ~{pairedFastq_1})
        fileBasename=${fileBasename%~{pairedSuffix_1}.*}

        # bwa alignment
        bwa mem ~{bwaParams} -t~{threadCount} asm.fa ~{pairedFastq_1} ~{pairedFastq_2} | samtools view -b -h > ${fileBasename}.bam
        samtools sort -@~{threadCount} -o ${fileBasename}.sorted.bam ${fileBasename}.bam
        du -s -BG ${fileBasename}.sorted.bam | sed 's/G.*//' > outputsize.txt
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
        File sortedBamFile = glob("*.sorted.bam")[0]
        Int fileSizeGB = read_int("outputsize.txt")
    }
}


