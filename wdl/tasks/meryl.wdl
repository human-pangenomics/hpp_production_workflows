version 1.0

import "extract_reads.wdl" as extractReads_t

workflow runMeryl {

    input {
        Array[File] childReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File? referenceFile
        Int threadCount
        String dockerImage
    }

    # actual work
    scatter (readFile in childReadsILM) {
        call extractReads_t.extractReads as childReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                threadCount=threadCount,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in maternalReadsILM) {
        call extractReads_t.extractReads as maternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                threadCount=threadCount,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsILM) {
        call extractReads_t.extractReads as paternalReadsExtracted {
            input:
                readFile=readFile,
                referenceFile=referenceFile,
                threadCount=threadCount,
                dockerImage=dockerImage
        }
    }
    call merylCount as childMerylCount {
        input:
            readFiles=childReadsExtracted.extractedRead,
            identifier="child",
            threadCount=threadCount,
            dockerImage=dockerImage
    }
    call merylCount as maternalMerylCount {
        input:
            readFiles=maternalReadsExtracted.extractedRead,
            identifier="maternal",
            threadCount=threadCount,
            dockerImage=dockerImage
    }
    call merylCount as paternalMerylCount {
        input:
            readFiles=paternalReadsExtracted.extractedRead,
            identifier="paternal",
            threadCount=threadCount,
            dockerImage=dockerImage
    }
    call merylHapmer as merylHapmer {
        input:
            childMerylDB=childMerylCount.merylDb,
            maternalMerylDB=maternalMerylCount.merylDb,
            paternalMerylDB=paternalMerylCount.merylDb,
            threadCount=threadCount,
            dockerImage=dockerImage
    }

	output {
		File childMerylDB = childMerylCount.merylDb
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
        Int threadCount
        String dockerImage
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
            meryl k=~{kmerSize} cpus=~{threadCount} count output reads$i.meryl $r
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
        cpu: threadCount
        memory: "8 GB"
        docker: dockerImage
    }
}


task merylHapmer {
    input {
        File maternalMerylDB
        File paternalMerylDB
        File childMerylDB
        Int threadCount
        String dockerImage
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
        tar xvf ~{childMerylDB} &
        wait

        # generate hapmers
        $MERQURY/trio/hapmers.sh maternal.meryl paternal.meryl child.meryl

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
        cpu: threadCount
        memory: "8 GB"
        docker: dockerImage
    }
}

