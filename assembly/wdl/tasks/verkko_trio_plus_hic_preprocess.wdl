version 1.0

import "../../../QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadsToGZ_t
import "../../../QC/wdl/tasks/meryl.wdl" as meryl_t
import "filter_hifi_adapter.wdl" as adapter_t

workflow verkko_preprocess_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Extracts/pull files for Verkko trio+hic run. Also creates hapmer DBs"
    }
    
    input {
        Array[File] input_hifi
        Array[File] input_nanopore
        Array[File] sample_illumina
        Array[File] maternal_illumina
        Array[File] paternal_illumina
        Array[File] hic

        
        File? referenceFasta        
        Boolean filterAdapters = true
        String excludeStringReadExtraction=""

        Int preemptible=1
        Int fileExtractionDiskSizeGB = 512
    }


    scatter (readFile in input_hifi) {
        call extractReadsToGZ_t.extractReadstoGZ as childReadsHiFiExtractedGz {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                excludeString=excludeStringReadExtraction,
                diskSizeGB=fileExtractionDiskSizeGB
        }

        if (filterAdapters){
            call adapter_t.cutadapt as filterAdapterHiFi {
                input:
                    readFastqGz = childReadsHiFiExtractedGz.extractedRead,
                    diskSizeGB = fileExtractionDiskSizeGB
            } 
        }
        File extracted_hifi = select_first([filterAdapterHiFi.filteredReadFastqGz, childReadsHiFiExtractedGz.extractedRead])
    }
       

    scatter (readFile in input_nanopore) {
        call extractReadsToGZ_t.extractReadstoGZ as extract_ultralong {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                excludeString=excludeStringReadExtraction,
                diskSizeGB=fileExtractionDiskSizeGB
        }
    }


    ## Create meryl hapmer DBs
    call meryl_t.runMeryl as meryl {
        input:
            compress         = true,
            sampleReadsILM   = sample_illumina,
            maternalReadsILM = maternal_illumina,
            paternalReadsILM = paternal_illumina,
            referenceFasta   = referenceFasta
    }

    Array[File] hic_localized = hic

    output {
        Array[File] ultralong_fq   = extract_ultralong.extractedRead
        Array[File] hifi_fq        = extracted_hifi
        File mat_compr_hapmer      = meryl.maternalHapmer
        File pat_compr_hapmer      = meryl.paternalHapmer
        File sample_compr_meryl_db = meryl.sampleMerylDB
        Array[File] out_hic        = hic_localized
    }
}

