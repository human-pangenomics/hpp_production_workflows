version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "sum.wdl" as sum_t

workflow runMeryl {

    input {
        Array[File] sampleReadsILM
        Array[File] maternalReadsILM
        Array[File] paternalReadsILM
        File? referenceFasta
        Int kmerSize = 21
        Int shardLinesPerFile = 256000000
        Int memSizeGB = 32
        Int threadCount = 32
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_merqury:latest"
    }

    # extract reads
    scatter (readFile in sampleReadsILM) {
        call extractReads_t.extractReads as sampleReadsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
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
                referenceFasta=referenceFasta,
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
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results (for union sum)
    call sum_t.sum as sampleReadSize {
        input:
            integers=sampleReadsExtracted.fileSizeGB
    }
    call sum_t.sum as maternalReadSize {
        input:
            integers=maternalReadsExtracted.fileSizeGB
    }
    call sum_t.sum as paternalReadSize {
        input:
            integers=paternalReadsExtracted.fileSizeGB
    }
    call sum_t.sum as allReadSize {
        input:
            integers=[sampleReadSize.value, maternalReadSize.value, paternalReadSize.value]
    }

    # shard reads
    scatter (readFile in sampleReadsExtracted.extractedRead) {
        call shardReads_t.shardReads as sampleShardReads {
            input:
                readFile=readFile,
                linesPerFile=shardLinesPerFile,
                threadCount=1,
                memSizeGB=4,
                diskSizeGB=fileExtractionDiskSizeGB*2,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in maternalReadsExtracted.extractedRead) {
        call shardReads_t.shardReads as maternalShardReads {
            input:
                readFile=readFile,
                linesPerFile=shardLinesPerFile,
                threadCount=1,
                memSizeGB=4,
                diskSizeGB=fileExtractionDiskSizeGB*2,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in paternalReadsExtracted.extractedRead) {
        call shardReads_t.shardReads as paternalShardReads {
            input:
                readFile=readFile,
                linesPerFile=shardLinesPerFile,
                threadCount=1,
                memSizeGB=4,
                diskSizeGB=fileExtractionDiskSizeGB*2,
                dockerImage=dockerImage
        }
    }

    # do the meryl counting
    scatter (readFile in flatten(sampleShardReads.shardedReads)) {
        call merylCount as sampleMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                threadCount=threadCount,
                memSizeGB=memSizeGB,
                diskSizeGB=sampleShardReads.fileSizeGB[0] * 4,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in flatten(maternalShardReads.shardedReads)) {
        call merylCount as maternalMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                threadCount=threadCount,
                memSizeGB=memSizeGB,
                diskSizeGB=maternalShardReads.fileSizeGB[0] * 4,
                dockerImage=dockerImage
        }
    }
    scatter (readFile in flatten(paternalShardReads.shardedReads)) {
        call merylCount as paternalMerylCount {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                threadCount=threadCount,
                memSizeGB=memSizeGB,
                diskSizeGB=paternalShardReads.fileSizeGB[0] * 4,
                dockerImage=dockerImage
        }
    }

    # do the meryl merging
    call merylUnionSum as sampleMerylUnionSum {
        input:
            merylCountFiles=sampleMerylCount.merylDb,
            identifier="sample",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=sampleReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylUnionSum as maternalMerylUnionSum {
        input:
            merylCountFiles=maternalMerylCount.merylDb,
            identifier="maternal",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=maternalReadSize.value * 4,
            dockerImage=dockerImage
    }
    call merylUnionSum as paternalMerylUnionSum {
        input:
            merylCountFiles=paternalMerylCount.merylDb,
            identifier="paternal",
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=paternalReadSize.value * 4,
            dockerImage=dockerImage
    }

    call merylHapmer as merylHapmer {
        input:
            sampleMerylDB=sampleMerylUnionSum.merylDb,
            maternalMerylDB=maternalMerylUnionSum.merylDb,
            paternalMerylDB=paternalMerylUnionSum.merylDb,
            threadCount=threadCount,
            memSizeGB=memSizeGB,
            diskSizeGB=allReadSize.value * 2,
            dockerImage=dockerImage
    }

	output {
		File sampleMerylDB = sampleMerylUnionSum.merylDb
		File maternalHapmer = merylHapmer.maternalHapmers
		File paternalHapmer = merylHapmer.paternalHapmers
		File hapmerImages = merylHapmer.hapmerImages
	}
}


task merylCount {
    input {
        File readFile
        Int kmerSize=21
        Int memSizeGB = 32
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
        ID=`basename ~{readFile} | sed 's/.gz$//' | sed 's/.f[aq]\(st[aq]\)*$//'`
        meryl k=~{kmerSize} threads=~{threadCount} memory=$((~{memSizeGB}-10)) count output $ID.meryl ~{readFile}

        # package
        tar cvf $ID.meryl.tar $ID.meryl

        # cleanup
        rm -rf $ID.meryl
	>>>
	output {
		File merylDb = glob("*.meryl.tar")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}


task merylUnionSum {
    input {
        Array[File] merylCountFiles
        String identifier
        Int memSizeGB = 64
        Int threadCount = 32
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
        mkdir extracted/
        cd extracted/
        for m in ~{sep=" " merylCountFiles} ; do
            tar xf $m &
        done
        wait
        cd ../

        # merge meryl dbs
        meryl union-sum output ~{identifier}.meryl extracted/*

        # package
        tar cvf ~{identifier}.meryl.tar ~{identifier}.meryl

        # cleanup
        rm -rf extracted/
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
        Int memSizeGB = 32
        Int threadCount = 32
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
        ln -s /bin/bash /usr/bin/sh

        # extract meryl dbs
        tar xvf ~{maternalMerylDB} &
        tar xvf ~{paternalMerylDB} &
        tar xvf ~{sampleMerylDB} &
        wait

        # generate hapmers
        bash $MERQURY/trio/hapmers.sh maternal.meryl paternal.meryl sample.meryl

        # package images
        tar czvf hapmers_img.tar.gz *.png *.hist

        # our desired destination is a softlink, need to move files to a real directory to tar
        mv maternal.hapmer.meryl maternal.tmp
        mkdir maternal.hapmers.meryl
        mv maternal.tmp/* maternal.hapmers.meryl/
        mv paternal.hapmer.meryl paternal.tmp
        mkdir paternal.hapmers.meryl
        mv paternal.tmp/* paternal.hapmers.meryl/
        tar cvf maternal.hapmers.meryl.tar maternal.hapmers.meryl &
        tar cvf paternal.hapmers.meryl.tar paternal.hapmers.meryl &
        wait
	>>>
	output {
		File maternalHapmers = "maternal.hapmers.meryl.tar"
		File paternalHapmers = "paternal.hapmers.meryl.tar"
		File hapmerImages = "hapmers_img.tar.gz"
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

