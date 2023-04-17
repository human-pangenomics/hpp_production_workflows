version 1.0

import "../../../QC/wdl/tasks/extract_reads.wdl" as extractReads_t
import "../../../QC/wdl/tasks/arithmetic.wdl" as arithmetic_t
import "filter_hifi_adapter.wdl" as adapter_t
import "filter_short_reads.wdl" as filter_short_reads_t

workflow runTrioHifiasm{
    input {
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        Array[File] childReadsONT=[]
        Int? homCov
        Int minOntReadLength=100000
        Int minHiFiReadLength=1
        String childID
        String? hifiasmExtraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
        Boolean filterAdapters
        Array[Float] offsetMem = [10, 0, 0]
        Array[Float] memCovRatios = [4.7, 3.8, 3.6]
        String excludeStringReadExtraction=""
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "quay.io/masri2019/hpp_hifiasm:0.18.5-r500"
        String zones = "us-west2-a"
    }

    scatter (readFile in childReadsHiFi) {
        call extractReads_t.extractReads as childReadsHiFiExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                excludeString=excludeStringReadExtraction,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
        # filter short HiFi reads
        call filter_short_reads_t.filterShortReads as extractLongHiFiReads{
            input:
                readFastq = childReadsHiFiExtracted.extractedRead,
                diskSizeGB = fileExtractionDiskSizeGB,
                minReadLength = minHiFiReadLength
        }
        if (filterAdapters){
            call adapter_t.cutadapt as filterAdapterHiFi {
                input:
                    readFastqGz = extractLongHiFiReads.longReadFastqGz,
                    diskSizeGB = fileExtractionDiskSizeGB
            } 
        }
        File hifiProcessedFastqGz = select_first([filterAdapterHiFi.filteredReadFastqGz, extractLongHiFiReads.longReadFastqGz])
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

    call arithmetic_t.sum as childReadHiFiSize {
        input:
            integers=extractLongHiFiReads.fileSizeGB
    }

    # if no ONT data is provided then it would be zero
    Int readULSize = select_first([childReadULSize.value, 0])

    call trioHifiasm as hifiasmStep1{
        input:
            childReadsHiFi=hifiProcessedFastqGz,
            homCov = homCov,
            childID=childID,
            extraOptions="--bin-only",
            #inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=ceil(memCovRatios[0] * homCov + offsetMem[0]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize.value + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    call trioHifiasm as hifiasmStep2{
        input:
            childReadsUL=extractUltraLongReads.longReadFastqGz, # optional argument
            paternalYak=paternalYak,
            maternalYak=maternalYak,
            homCov = homCov,
            childID=childID,
            extraOptions="--bin-only",
            inputBinFilesTarGz=hifiasmStep1.outputBinFiles,
            memSizeGB=ceil(memCovRatios[1] * homCov + offsetMem[1]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize.value + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    call trioHifiasm as hifiasmStep3{
        input:
            childReadsUL=extractUltraLongReads.longReadFastqGz, # optional argument
            homCov = homCov,
            childID=childID,
            extraOptions="",
            inputBinFilesTarGz=hifiasmStep2.outputBinFiles,
            memSizeGB=ceil(memCovRatios[2] * homCov + offsetMem[2]),
            threadCount=threadCount,
            diskSizeGB= floor((childReadHiFiSize.value + readULSize) * 2.5) + 64,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones = zones
    }
    
    output {
        File outputPaternalGfa = hifiasmStep3.outputPaternalGfa
        File outputMaternalGfa = hifiasmStep3.outputMaternalGfa
        File outputPaternalContigGfa = hifiasmStep3.outputPaternalContigGfa
        File outputMaternalContigGfa = hifiasmStep3.outputMaternalContigGfa
        File outputRawUnitigGfa = hifiasmStep3.outputRawUnitigGfa
        File outputBinFiles = hifiasmStep3.outputBinFiles
    }
}

# hifiasm steps
# 1st: pass HiFi and both yak files / not need to pass UL (extraOptions="--bin-only")
# 2nd: pass UL and both yak files (extraOptions="--bin-only")
# 3rd: pass UL no need to pass yak files 
task trioHifiasm {
    input{
        File? paternalYak
        File? maternalYak
        Array[File]? childReadsHiFi
        Array[File]? childReadsUL
        Int? homCov
        String childID
	String? extraOptions
        File? inputBinFilesTarGz
        File? referenceFasta
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

        ## run trio hifiasm https://github.com/chhylp123/hifiasm
        # If ONT ultra long reads are provided
        if [[ -n "~{sep="" childReadsUL}" ]]; then
            if [[ -n "~{paternalYak}" ]]; then 
                # Run step 2
                hifiasm ~{extraOptions} -o ~{childID} --ul ~{sep="," childReadsUL} --hom-cov ~{homCov} -t~{threadCount} -1 ~{paternalYak} -2 ~{maternalYak} fake.fq
            else
                # Keep only necessary bin files for step 3
                mkdir kept_bin_files
                mv *.ec.bin *.ovlp.reverse.bin *.ovlp.source.bin *.hap1.phase.bin *.hap2.phase.bin *.ul.ovlp.bin kept_bin_files
                rm -rf *.bin
                mv kept_bin_files/* .
                rm -rf kept_bin_files 
                
                # Run step 3
                hifiasm ~{extraOptions} -o ~{childID} --ul ~{sep="," childReadsUL} --hom-cov ~{homCov} -t~{threadCount} -3 ~{childID}.hap1.phase.bin -4 ~{childID}.hap2.phase.bin fake.fq
            fi
        else  
            # Run step 1
            hifiasm ~{extraOptions} -o ~{childID} -t~{threadCount} ~{sep=" " childReadsHiFi}
        fi

        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.pat.contig_gfa
        mkdir ~{childID}.mat.contig_gfa
        mkdir ~{childID}.binFiles

        # before hardlinking gfa files to the corresponding directory make sure they exist
        # first and second step of hifiasm does not output gfa files
        if [[ -n $(find . -maxdepth 1 -iname "*.hap1.p_ctg.*") ]];
        then
            ln ~{childID}.dip.r_utg.* ~{childID}.raw_unitig_gfa
            ln *.hap1.p_ctg.* ~{childID}.pat.contig_gfa
            ln *.hap2.p_ctg.* ~{childID}.mat.contig_gfa
        else
            # To avoid making a separete new task for steps 1 and 2 of hifiasm
            # we make empty gfa files since we cannot have optional outputs in wdl
            touch empty.hap1.p_ctg.gfa
            touch empty.hap2.p_ctg.gfa
        fi

        ln *.bin ~{childID}.binFiles
        
        
        # make archives
        tar -cf ~{childID}.raw_unitig_gfa.tar ~{childID}.raw_unitig_gfa
        tar -cf ~{childID}.pat.contig_gfa.tar ~{childID}.pat.contig_gfa
        tar -cf ~{childID}.mat.contig_gfa.tar ~{childID}.mat.contig_gfa
        tar -cf ~{childID}.binFiles.tar ~{childID}.binFiles
        
        # compress
        pigz -p~{threadCount} ~{childID}.raw_unitig_gfa.tar
        pigz -p~{threadCount} ~{childID}.pat.contig_gfa.tar
        pigz -p~{threadCount} ~{childID}.mat.contig_gfa.tar
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
        File outputPaternalGfa = glob("*.hap1.p_ctg.gfa")[0]
        File outputMaternalGfa = glob("*.hap2.p_ctg.gfa")[0]
        File outputPaternalContigGfa = "~{childID}.pat.contig_gfa.tar.gz"
        File outputMaternalContigGfa = "~{childID}.mat.contig_gfa.tar.gz"
        File outputRawUnitigGfa = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFiles = "~{childID}.binFiles.tar.gz"
    }
}

