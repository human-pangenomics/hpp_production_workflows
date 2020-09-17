version 1.0

import "extract_reads.wdl" as extractReads_t

workflow runMeryl {

    input {
        Array[File] sampleReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File? referenceFile
        Int memSizeGB = 64
        Int threadCount = 16
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_merqury:latest"
    }

    # extract reads
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in maternalReadsILM) {
        call extractReads_t.extractReads as maternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsILM) {
        call extractReads_t.extractReads as paternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results
    call sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }
    call sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call sum as allReadSize {
        input:
            integers=[sampleReadSize.value, maternalReadSize.value, paternalReadSize.value]
    }

    # do the actual meryl work
    call merylCount as sampleMerylCount {
        input:
            readFiles=sampleReadsExtracted.extractedRead,
            identifier="sample",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=sampleReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylCount as maternalMerylCount {
        input:
            readFiles=maternalReadsExtracted.extractedRead,
            identifier="maternal",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=maternalReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylCount as paternalMerylCount {
        input:
            readFiles=paternalReadsExtracted.extractedRead,
            identifier="paternal",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=maternalReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylHapmer as merylHapmer {
        input:
            sampleMerylDB=sampleMerylCount.merylDb,
            maternalMerylDB=maternalMerylCount.merylDb,
            paternalMerylDB=paternalMerylCount.merylDb,
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=allReadSize.value,
            dockerImage=dockerImage
    }

	output {
		File sampleMerylDB = sampleMerylCount.merylDb
		File maternalHapmer = merylHapmer.maternalHapmers
		File paternalHapmer = merylHapmer.paternalHapmers
		File hapmerImage = merylHapmer.hapmerImage
	}
}


task merylCount {
    input {
        Array[File] readFiles
        String identifier
        Int kmerSize=21
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # generate meryl db for each read
        i=0
        for r in ~{sep=" " readFiles} ; do
            meryl k=~{kmerSize} threads=~{threadCount} memory=$((~{memSizeGB}-12)) count output reads$i.meryl $r
            i=$(($i + 1))
        done

        # merge meryl dbs
        meryl cpus=~{threadCount} union-sum output ~{identifier}.meryl reads*.meryl

        # package
        tar vf ~{identifier}.meryl.tar ~{identifier}.meryl
	>>>
	output {
		File merylDb = identifier + ".meryl.tar"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}


task merylHapmer {
    input {
        File maternalMerylDB
        File paternalMerylDB
        File sampleMerylDB
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_merqury:latest"
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
        OMP_NUM_THREADS=~{threadCount}

        # extract meryl dbs
        tar xvf ~{maternalMerylDB} &
        tar xvf ~{paternalMerylDB} &
        tar xvf ~{sampleMerylDB} &
        wait

        # generate hapmers
        $MERQURY/trio/hapmers.sh maternal.meryl paternal.meryl sample.meryl

        # package
        tar vf mat.hapmers.meryl.tar mat.hapmers.meryl &
        tar vf pat.hapmers.meryl.tar pat.hapmers.meryl &
        wait
	>>>
	output {
		File maternalHapmers = "mat.hapmers.meryl.tar"
		File paternalHapmers = "pat.hapmers.meryl.tar"
		File hapmerImage = "inherited_hapmers.png"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}



task sum {
    input {
        Array[Int] integers
    }

    command <<<
        echo $((~{sep="+" integers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        docker: "tpesout/hpp_base:latest"
    }
}

