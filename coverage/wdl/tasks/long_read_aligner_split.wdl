version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "merge_bams.wdl" as mergeBams_t
import "read_set_splitter.wdl" as readSetSplitter_t
import "long_read_aligner.wdl" as longReadAligner_t

workflow longReadAlignmentSplit {
    input {
        String aligner="winnowmap"
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        Int splitNumber = 16
        File assembly
        File? referenceFasta
        Int preemptible=2
        Int extractReadsDiskSize=256
        String zones
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="tpesout/hpp_base:latest"
        }
    }
    call readSetSplitter_t.readSetSplitter {
        input:
            readFastqs = extractReads.extractedRead,
            splitNumber = splitNumber
    }
    call arithmetic_t.sum as readSize {
        input:
            integers=extractReads.fileSizeGB
    }
    scatter (readFastq in readSetSplitter.splitReadFastqs) {
         ## align reads to the assembly
         call longReadAligner_t.alignment{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=assembly,
                 readFastq_or_queryAssembly = readFastq,
                 diskSize = floor(readSize.value / splitNumber) * 4,
                 preemptible = preemptible,
                 zones = zones
        }
    }
    call arithmetic_t.sum as bamSize {
        input:
            integers=alignment.fileSizeGB
    }

    ## merge the bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = "${sampleName}.${sampleSuffix}",
            sortedBamFiles = alignment.sortedBamFile,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5),
            preemptible = preemptible,
            zones = zones
    }
    output {
        File sortedBamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }
}

