version 1.0

import "extract_reads.wdl" as extractReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runTrioHifiasm{
    input {
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        String childID
        File? inputBinFilesTarGz
        File? referenceFasta
        Int memSizeGB
        Int threadCount
        Int preemptible
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "quay.io/masri2019/hpp_hifiasm:latest"
        String zones
    }

    scatter (readFile in childReadsHiFi) {
        call extractReads_t.extractReads as childReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    call arithmetic_t.sum as childReadSize {
        input:
            integers=childReadsExtracted.fileSizeGB
    }

    call trioHifiasm {
        input:
            paternalYak=paternalYak,
            maternalYak=maternalYak,
            childReadsHiFi=childReadsHiFi,
            childID=childID,
            inputBinFilesTarGz=inputBinFilesTarGz,
            memSizeGB=memSizeGB,
            threadCount=threadCount,
            diskSizeGB=childReadSize.value * 2,
            preemptible=preemptible,
            dockerImage=dockerImage,
            zones=zones
    }
    output {
        File outputPaternalGfa = trioHifiasm.outputPaternalGfa
        File outputMaternalGfa = trioHifiasm.outputMaternalGfa
        File outputPaternalContigGfaTarGz = trioHifiasm.outputPaternalContigGfaTarGz
        File outputMaternalContigGfaTarGz = trioHifiasm.outputMaternalContigGfaTarGz
        File outputRawUnitigGfaTarGz = trioHifiasm.outputRawUnitigGfaTarGz
        File outputBinFilesTarGz = trioHifiasm.outputBinFilesTarGz
    }
}

task trioHifiasm {
    input{
        File paternalYak
        File maternalYak
        Array[File] childReadsHiFi
        String childID
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

        ## if bin files are given we have to extract them in the directory where hifiasm is being run
        ## this enables hifiasm to skip the time-consuming process of finding read overlaps
        if [ ! -v ~{inputBinFilesTarGz} ]; then
            tar -xzf ~{inputBinFilesTarGz} --strip-components 1
            rm -rf ~{inputBinFilesTarGz}
        fi

        ## it is not possible to pipe multiple fastq files to hifiasm
        cat ~{sep=" " childReadsHiFi} > ~{childID}.fastq
        rm -rf ~{sep=" " childReadsHiFi}

        ## run trio hifiasm https://github.com/chhylp123/hifiasm
        hifiasm -o ~{childID} -t~{threadCount} -1 ~{paternalYak} -2 ~{maternalYak} ~{childID}.fastq
        
        #Move bin and gfa files to saparate folders and compress them 
        mkdir ~{childID}.raw_unitig_gfa
        mkdir ~{childID}.pat.contig_gfa
        mkdir ~{childID}.mat.contig_gfa
        mkdir ~{childID}.binFiles
        
        ln ~{childID}.dip.r_utg.* ~{childID}.raw_unitig_gfa
        ln ~{childID}.hap1.p_ctg.* ~{childID}.pat.contig_gfa
        ln ~{childID}.hap2.p_ctg.* ~{childID}.mat.contig_gfa
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
        zones : zones
    }

    output {
        File outputPaternalGfa = "~{childID}.hap1.p_ctg.gfa"
        File outputMaternalGfa = "~{childID}.hap2.p_ctg.gfa"
        File outputPaternalContigGfaTarGz = "~{childID}.pat.contig_gfa.tar.gz"
        File outputMaternalContigGfaTarGz = "~{childID}.mat.contig_gfa.tar.gz"
        File outputRawUnitigGfaTarGz = "~{childID}.raw_unitig_gfa.tar.gz"
        File outputBinFilesTarGz = "~{childID}.binFiles.tar.gz"
    }
}

