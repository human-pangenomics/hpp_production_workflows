version 1.0

import "extract_reads.wdl" as extractReads_t
import "arithmetic.wdl" as arithmetic_t
import "merge_bams.wdl" as mergeBams_t
import "read_set_splitter.wdl" as readSetSplitter_t
import "long_read_aligner.wdl" as longReadAligner_t

workflow longReadAlignmentScattered {
    input {
        # aligner:
        #     Can be either minimap2, winnowmap and veritymap
        String aligner="winnowmap"
        # preset:
        #     For winnowmap-v2.03 -> map-pb/map-ont/asm5/asm10
        #     For minimap2-v2.24 -> map-pb/map-hifi/map-ont/asm5/asm10
        #     For veritymap v2.1.2-alpha -> hifi-haploid/hifi-haploid-complete/hifi-diploid/ont-haploid-complete
        String preset
        String sampleName
        String sampleSuffix
        Array[File] readFiles
        Int splitNumber = 4
        File assembly
        File? referenceFasta
        # options:
        #     For winnowmap/minimap2 recommended "--eqx --cs -Y -L" (if diploid assembly "-I8g")
        #     For veritymap "--careful" can be used but not recommended for whole-genome assembly
        #     Only winnowmap/minimap2 supports methylation tags in reads (parameter should contain "-TMl,Mm")
        String options=""
        # fastqOptions:
        #     For reads with methylation "-y"
        String fastqOptions=""
        Int kmerSize = 15
        Int preemptible=2
        Int extractReadsDiskSize=512
        String zones="us-west2-a"
        String dockerImage="mobinasri/long_read_aligner:v0.2"
    }

    scatter (readFile in readFiles) {
        call extractReads_t.extractReads as extractReads {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=extractReadsDiskSize,
                dockerImage="mobinasri/bio_base:v0.2"
        }
    }

   call arithmetic_t.sum as readSize {
        input:
            integers=extractReads.fileSizeGB
    }

    call readSetSplitter_t.readSetSplitter {
        input:
            readFastqs = extractReads.extractedRead,
            splitNumber = splitNumber,
            diskSize = floor(readSize.value * 2.5)
    }


   scatter (readFastqAndSize in zip(readSetSplitter.splitReadFastqs, readSetSplitter.readSizes)) {
         ## align reads to the assembly
         call longReadAligner_t.alignmentBam as alignment{
             input:
                 aligner =  aligner,
                 preset = preset,
                 refAssembly=assembly,
                 readFastq_or_queryAssembly = readFastqAndSize.left,
                 diskSize = 64 + floor(readFastqAndSize.right) * 6,
                 preemptible = preemptible,
                 zones = zones,
                 options = options,
                 kmerSize = kmerSize,
                 dockerImage=dockerImage
         }
    }

    call arithmetic_t.sum as bamSize {
        input:
            integers = alignment.fileSizeGB
    }


    ## merge the bam files
    call mergeBams_t.merge as mergeBams{
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            sortedBamFiles = alignment.sortedBamFile,
            # runtime configurations
            diskSize = floor(bamSize.value * 2.5) + 32,
            preemptible = preemptible,
            zones = zones,
            dockerImage=dockerImage
    }

    output {
        File bamFile = mergeBams.mergedBam
        File baiFile = mergeBams.mergedBai
    }
}
