version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "filter_hifi_adapter.wdl" as adapter_t
import "filter_short_reads.wdl" as filter_short_reads_t
import "hifiasm_cutadapt_hic_multistep.wdl" as hifiasm_hic_t


workflow runHiCHifiasmStep2{
    input {
        File inputBinFilesTarGz # it should be provided from a previous run (either trio- or hic-based)
        Array[File] childReadsONT
        Int homCov
        Int minOntReadLength=100000
        String childID
        String? hifiasmExtraOptions
        File? referenceFasta
        Boolean filterAdapters
        Array[Float] offsetMem = [10, 0, 0]
        Array[Float] memCovRatios = [4.7, 3.8, 3.6]
        String excludeStringReadExtraction=""
	File fakeFastq = "gs://masri/hprc/fake.fq" 
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "quay.io/masri2019/hpp_hifiasm:0.18.5-r500"
        String zones = "us-west2-a"
    }

    # if ONT reads are provided
    if (length(childReadsONT) != 0){
        scatter (readFile in childReadsONT) {
            call extractReads_t.extractReads as childReadsOntExtracted {
                input:
                    readFile=readFile,
                    referenceFasta=referenceFasta,
                    memSizeGB=4,
                    threadCount=4,
                    diskSizeGB=fileExtractionDiskSizeGB,
                    dockerImage=dockerImage
            }
            # filter ONT reads to get UL reads
            call filter_short_reads_t.filterShortReads as extractUltraLongReads{
                input:
                    readFastq = childReadsOntExtracted.extractedRead,
                    diskSizeGB = fileExtractionDiskSizeGB,
                    minReadLength = minOntReadLength
            }
         }
         call arithmetic_t.sum as childReadULSize {
             input:
                 integers=extractUltraLongReads.fileSizeGB
         }
    }

    # if no ONT data is provided then it would be zero
    Int readULSize = select_first([childReadULSize.value, 0])

    call hifiasm_hic_t.hicHifiasm as hifiasmStep2{
        input:
            childReadsHiFi=[fakeFastq],
            childReadsUL=extractUltraLongReads.longReadFastqGz,
            homCov = homCov,
            childID=childID,
            extraOptions="--bin-only",
            inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=ceil(memCovRatios[1] * homCov + offsetMem[1]),
            threadCount=threadCount,
            diskSizeGB= floor(readULSize * 2.5) + 512,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    
    output {
        File outputBinFiles = hifiasmStep2.outputBinFiles
    }
}

