version 1.0

import "../../../QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadsToGZ_t
import "filter_hifi_adapter.wdl" as adapter_t
import "seqkit_filter_fastq.wdl" as seqkit_filter_wf

workflow runHiCHifiasm{
    input {
        Array[File] childReadsHiFi
        Array[File] childReadsONT=[]
        Array[File] childReadsHiC1
        Array[File] childReadsHiC2        
        Int? homCov
        Int minOntReadLength=50000
        String childID
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        Boolean filterAdapters
        Int? min_ont_qscore        
        Array[Float] offsetMem = [40, 20, 20]
        Array[Float] memCovRatios = [4.7, 3.8, 3.6]
        String excludeStringReadExtraction=""
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "humanpangenomics/hifiasm@sha256:50483d0757e5ce1eaf68f33da8ac8087f8b40f1d1310f60100945fc3d56721c5"
        String zones = "us-west2-a"
    }

    scatter (readFile in childReadsHiFi) {
        call extractReadsToGZ_t.extractReadstoGZ as childReadsHiFiExtractedGz {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                excludeString=excludeStringReadExtraction,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
        if (filterAdapters){
            call adapter_t.cutadapt as filterAdapterHiFi {
                input:
                    readFastqGz = childReadsHiFiExtractedGz.extractedRead,
                    diskSizeGB = fileExtractionDiskSizeGB
            } 
        }
        File hifiProcessedFastqGz = select_first([filterAdapterHiFi.filteredReadFastqGz, childReadsHiFiExtractedGz.extractedRead])
    }

    Float childReadHiFiSize = size(hifiProcessedFastqGz, 'GB')


    # if ONT reads are provided
    if (length(childReadsONT) != 0){
        scatter (readFile in childReadsONT) {
            call extractReadsToGZ_t.extractReadstoGZ as childReadsOntExtractedGz {
                input:
                    readFile=readFile,
                    referenceFasta=referenceFasta,
                    memSizeGB=4,
                    threadCount=4,
                    excludeString=excludeStringReadExtraction,
                    diskSizeGB=fileExtractionDiskSizeGB,
                    dockerImage=dockerImage
            }
            
            # filter ONT reads for quality
            if (defined(min_ont_qscore)) {
                call seqkit_filter_wf.filter_fastq as qscore_filter_ont {
                    input:
                        input_fastq  = childReadsOntExtractedGz.extractedRead,
                        min_size     = minOntReadLength,
                        min_q        = min_ont_qscore
                }
            }
            File ont_reads = select_first([qscore_filter_ont.filteredFastq, childReadsOntExtractedGz.extractedRead])
         }
         ## get size of all ultralong reads
         Float ul_size_all = size(ont_reads, 'GB')
    }
    # if no ONT data is provided then it would be zero
    Float readULSize = select_first([ul_size_all, 0])    


    call hicHifiasm as hifiasmStep1{
        input:
            childReadsHiFi=hifiProcessedFastqGz,
            homCov = homCov,
            childID=childID,
            extraOptions="--bin-only --telo-m CCCTAA",
            memSizeGB=ceil(memCovRatios[0] * select_first([homCov, 0]) + offsetMem[0]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    call hicHifiasm as hifiasmStep2{
        input:
            childReadsUL=ont_reads, 
            homCov = homCov,
            minOntReadLength = minOntReadLength,
            childID=childID,
            extraOptions="--bin-only --telo-m CCCTAA",
            inputBinFilesTarGz=hifiasmStep1.outputBinFiles,
            memSizeGB=ceil(memCovRatios[1] * select_first([homCov, 0]) + offsetMem[1]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 128,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    call hicHifiasm as hifiasmStep3{
        input:
            childReadsUL=ont_reads,
            childReadsHiC1=childReadsHiC1,
            childReadsHiC2=childReadsHiC2,
            homCov = homCov,
            minOntReadLength = minOntReadLength,
            childID=childID,
            extraOptions="--telo-m CCCTAA",
            inputBinFilesTarGz=hifiasmStep2.outputBinFiles,
            memSizeGB=ceil(memCovRatios[2] * select_first([homCov, 0]) + offsetMem[2]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize + readULSize) * 2.5) + 128,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    
    output {
        File outputHaplotype1Gfa = hifiasmStep3.outputHaplotype1Gfa
        File outputHaplotype2Gfa = hifiasmStep3.outputHaplotype2Gfa
        File outputHaplotype1ContigGfa = hifiasmStep3.outputHaplotype1ContigGfa
        File outputHaplotype2ContigGfa = hifiasmStep3.outputHaplotype2ContigGfa
        File outputRawUnitigGfa = hifiasmStep3.outputRawUnitigGfa
        File outputBinFiles = hifiasmStep3.outputBinFiles
    }
}

# hifiasm HiC steps
# 1st: pass HiFi; no need to pass UL and HiC (extraOptions="--bin-only")
# 2nd: pass UL and fake HiFi, not need to pass HiC (extraOptions="--bin-only")
# 3rd: pass UL, HiC and fake HiFi
task hicHifiasm {
    input{
        Array[File]? childReadsHiC1
        Array[File]? childReadsHiC2
        Array[File]? childReadsHiFi
        Array[File]? childReadsUL
        Int? homCov
        Int minOntReadLength = 50000
        String childID
        String? extraOptions
        File? inputBinFilesTarGz
        # runtime configurations
        Int memSizeGB
        Int threadCount
        Int diskSizeGB
        Int preemptible
        String dockerImage
        String zones
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## if bin files are given we have to extract them in the directory where hifiasm is being run
        ## this enables hifiasm to skip the time-consuming process of finding read overlaps
        if [ ! -v ~{inputBinFilesTarGz} ]; then
            tar -xzf ~{inputBinFilesTarGz} --strip-components 1
            rm -rf ~{inputBinFilesTarGz} || true
        fi

        # make a fake fastq to use for step 2 and 3
        printf "@fake\nA\n+\nI\n" > fake.fq

        ## run hic hifiasm https://github.com/chhylp123/hifiasm
        # If ONT ultra long reads are provided
        if [[ -n "~{sep="" childReadsUL}" ]]; then
            if [[ -n "~{sep="" childReadsHiC1}" ]]; then                 
                # Keep only necessary bin files for step 3
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin *.ul.ovlp.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files

                # Run step 3
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    --ul ~{sep="," childReadsUL} \
                    --hom-cov ~{homCov} \
                    --ul-cut ~{minOntReadLength} \
                    --dual-scaf \
                    -t~{threadCount} \
                    --h1 "~{sep="," childReadsHiC1}" \
                    --h2 "~{sep="," childReadsHiC2}"  \
                    fake.fq
            else
                # Keep only necessary bin files for step 2
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files

                # Run step 2
                hifiasm \
                    ~{extraOptions} \
                    -o ~{childID} \
                    --ul ~{sep="," childReadsUL} \
                    --ul-cut ~{minOntReadLength} \
                    --hom-cov ~{homCov} \
                    -t~{threadCount} \
                    fake.fq
            fi
        else  
            # Run step 1
            hifiasm \
                ~{extraOptions} \
                -o ~{childID} \
                -t~{threadCount} \
                ~{sep=" " childReadsHiFi}
        fi

        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.hap1.contig_gfa
        mkdir ~{childID}.hap2.contig_gfa
        mkdir ~{childID}.binFiles

        # before hardlinking gfa files to the corresponding directory make sure they exist
        # first and second step of hifiasm does not output gfa files
        if [[ -n $(find . -maxdepth 1 -iname "*.hap1.p_ctg.*") ]];
        then
            ln *.r_utg.* ~{childID}.raw_unitig_gfa
            ln *.hap1.p_ctg.* ~{childID}.hap1.contig_gfa
            ln *.hap2.p_ctg.* ~{childID}.hap2.contig_gfa
        else
            # To avoid making a separete new task for steps 1 and 2 of hifiasm
            # we make empty gfa files since we cannot have optional outputs in wdl
            touch empty.hap1.p_ctg.gfa
            touch empty.hap2.p_ctg.gfa
        fi

        ln *.bin ~{childID}.binFiles
        
        
        # make archives
        tar -cf ~{childID}.raw_unitig_gfa.tar ~{childID}.raw_unitig_gfa
        tar -cf ~{childID}.hap1.contig_gfa.tar ~{childID}.hap1.contig_gfa
        tar -cf ~{childID}.hap2.contig_gfa.tar ~{childID}.hap2.contig_gfa
        tar -cf ~{childID}.binFiles.tar ~{childID}.binFiles
        
        # compress
        pigz -p~{threadCount} ~{childID}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{childID}.hap1.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.hap2.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.binFiles.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        cpuPlatform: "Intel Cascade Lake"
        zones : zones
    }

    output {
        File outputHaplotype1Gfa = glob("*.hap1.p_ctg.gfa")[0]
        File outputHaplotype2Gfa = glob("*.hap2.p_ctg.gfa")[0]
        File outputHaplotype1ContigGfa = "~{childID}.hap1.contig_gfa.tar.gz"
        File outputHaplotype2ContigGfa = "~{childID}.hap2.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{childID}.binFiles.tar.gz"
    }
}

