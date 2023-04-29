version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "filter_hifi_adapter.wdl" as adapter_t
import "filter_short_reads.wdl" as filter_short_reads_t
import "hifiasm_cutadapt_hic_multistep.wdl" as hifiasm_hic_t
import "gfatools.wdl" as gfatools_t


workflow runHiCHifiasmStep3{
    input {
        File inputBinFilesTarGz # it should be provided from a previous run (hic step2)
        Array[File] childReadsONT
        Array[File] childReadsHiC1
        Array[File] childReadsHiC2        
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
                    dockerImage="mobinasri/bio_base:v0.2"
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

    call hifiasm_hic_t.hicHifiasm as hifiasmStep3{
        input:
            childReadsHiFi = [fakeFastq],
            childReadsUL = extractUltraLongReads.longReadFastqGz,
            childReadsHiC1 = childReadsHiC1,
            childReadsHiC2 = childReadsHiC2,
            homCov = homCov,
            childID=childID,
            extraOptions="",
            inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=ceil(memCovRatios[1] * homCov + offsetMem[1]),
            threadCount=threadCount,
            diskSizeGB= floor(readULSize * 2.5) + 1024,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }

    call gfatools_t.phasedGFAs2Fasta as gfa2fasta{
        input:
            paternalGfa = hifiasmStep3.outputHaplotype1Gfa,
            maternalGfa = hifiasmStep3.outputHaplotype2Gfa,
            patSuffix = "hap1",
            matSuffix = "hap2",
            childID = childID
    }
 
    output {
        File haplotype1FastaGz = gfa2fasta.outputPaternalFastaGz
        File haplotype2FastaGz = gfa2fasta.outputMaternalFastaGz
        File outputHaplotype1ContigGfa = hifiasmStep3.outputHaplotype1ContigGfa
        File outputHaplotype2ContigGfa = hifiasmStep3.outputHaplotype2ContigGfa
        File outputRawUnitigGfa = hifiasmStep3.outputRawUnitigGfa
        File outputBinFiles = hifiasmStep3.outputBinFiles
    }
}

