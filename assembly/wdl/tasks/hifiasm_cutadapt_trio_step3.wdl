version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "filter_hifi_adapter.wdl" as adapter_t
import "filter_short_reads.wdl" as filter_short_reads_t
import "hifiasm_cutadapt_multistep.wdl" as hifiasm_trio_t
import "gfatools.wdl" as gfatools_t


workflow runTrioHifiasmStep3{
    input {
        File inputBinFilesTarGz # it should be provided from a previous run (trio step2)
        Array[File] childReadsONT
        File paternalYak
        File maternalYak        
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
        String dockerImage = "quay.io/masri2019/hpp_hifiasm:0.19.0"
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

    call hifiasm_trio_t.trioHifiasm as hifiasmStep3{
        input:
            paternalYak=paternalYak,
            maternalYak=maternalYak,
            childReadsHiFi=[fakeFastq],
            childReadsUL=extractUltraLongReads.longReadFastqGz, # optional argument
            homCov = homCov,
            childID=childID,
            extraOptions="",
            inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=ceil(memCovRatios[2] * homCov + offsetMem[2]),
            threadCount=threadCount,
            diskSizeGB= floor((readULSize) * 2.5) + 1024,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }


    call gfatools_t.phasedGFAs2Fasta as gfa2fasta{
        input:
            paternalGfa = hifiasmStep3.outputPaternalGfa,
            maternalGfa = hifiasmStep3.outputMaternalGfa,
            childID = childID
    }
 
    output {
        File paternalFastaGz = gfa2fasta.outputPaternalFastaGz
        File maternalFastaGz = gfa2fasta.outputMaternalFastaGz
        File outputPaternalContigGfa = hifiasmStep3.outputPaternalContigGfa
        File outputMaternalContigGfa = hifiasmStep3.outputMaternalContigGfa
        File outputRawUnitigGfa = hifiasmStep3.outputRawUnitigGfa
        File outputBinFiles = hifiasmStep3.outputBinFiles
    }
}

