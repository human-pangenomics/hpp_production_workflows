version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runMash {

    input {
        Array[File] reads
        File? referenceFasta
        Int kmerSize = 21
        Int shardLinesPerFile = 256000000
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_mash:latest"
    }

    # extract reads
    scatter (readFile in reads) {
        call extractReads_t.extractReads as readsExtracted {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }
    }

    # get file size of results
    call arithmetic_t.max as maxReadSize {
        input:
            integers=readsExtracted.fileSizeGB
    }

#    # shard reads
#    scatter (readFile in sampleReadsExtracted.extractedRead) {
#        call shardReads_t.shardReads as sampleShardReads {
#            input:
#                readFile=readFile,
#                linesPerFile=shardLinesPerFile,
#                threadCount=1,
#                memSizeGB=4,
#                diskSizeGB=fileExtractionDiskSizeGB*2,
#                dockerImage=dockerImage
#        }
#    }
#    scatter (readFile in maternalReadsExtracted.extractedRead) {
#        call shardReads_t.shardReads as maternalShardReads {
#            input:
#                readFile=readFile,
#                linesPerFile=shardLinesPerFile,
#                threadCount=1,
#                memSizeGB=4,
#                diskSizeGB=fileExtractionDiskSizeGB*2,
#                dockerImage=dockerImage
#        }
#    }
#    scatter (readFile in paternalReadsExtracted.extractedRead) {
#        call shardReads_t.shardReads as paternalShardReads {
#            input:
#                readFile=readFile,
#                linesPerFile=shardLinesPerFile,
#                threadCount=1,
#                memSizeGB=4,
#                diskSizeGB=fileExtractionDiskSizeGB*2,
#                dockerImage=dockerImage
#        }
#    }
#    scatter (readFile in flatten(maternalShardReads.shardedReads)) {
#        call merylCount as maternalMerylCount {
#            input:
#                readFile=readFile,
#                kmerSize=kmerSize,
#                threadCount=threadCount,
#                memSizeGB=memSizeGB,
#                diskSizeGB=maternalShardReads.fileSizeGB[0] * 4,
#                dockerImage=dockerImage
#        }
#    }

    # sketch
    scatter (readFile in readsExtracted.extractedRead) {
        call mashSketch {
            input:
                readFile=readFile,
                kmerSize=kmerSize,
                diskSizeGB=maxReadSize.value + 10,
                dockerImage=dockerImage
        }
    }

    # merge the paste
    call mashPaste  {
        input:
            sketches=mashSketch.sketch,
            identifier="all",
            dockerImage=dockerImage
    }

    # final results
    call mashDistPlot as allMashDistPlot {
        input:
            querySketch=mashPaste.paste,
            referenceSketch=mashPaste.paste,
            dockerImage=dockerImage
    }

    scatter (readFile in readsExtracted.extractedRead) {
        call mashScreen as allMashScreen {
            input:
                readFile=readFile,
                diskSizeGB=maxReadSize.value + 10,
                dockerImage=dockerImage
        }
    }

	output {
		File distTable = allMashDistPlot.table
		File distPlot = allMashDistPlot.plot
		Array[File] screenResult = allMashScreen.screenOut
	}
}


task mashSketch {
    input {
        File readFile
        Int kmerSize=21
        Int sketchSize=10000
        String extraArgs="-r -m 1"
        Int memSizeGB = 2
        Int threadCount = 2
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_mash:latest"
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

        # get output name
        OUTPUT=$(basename ~{readFile}).msh

        # sketch
        mash sketch -k ~{kmerSize} -s ~{sketchSize} ~{extraArgs} -o $OUTPUT ~{readFile}

	>>>
	output {
		File sketch = glob("*.msh")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}


task mashPaste {
    input {
        Array[File] sketches
        String identifier
        Int memSizeGB = 2
        Int threadCount = 2
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_mash:latest"
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

        # ensure identifier ends with .msh
        ID=`echo ~{identifier} | sed 's/.msh$//' | sed 's/$/.msh/'`

        # paste together
        mash paste $ID ~{sep=" " sketches}

	>>>
	output {
		File paste = glob("*.msh")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}


task mashDistPlot {
    input {
        File querySketch
        File referenceSketch
        String extraArguments=""
        Int memSizeGB = 8
        Int threadCount = 4
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_mash:latest"
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

        # this is a hack!
        mv /root/bin/R_3.6.3/bin/Rscript /root/bin/R_3.6.3/bin/.Rscript && \
        echo -e '#!/bin/bash\nxvfb-run -s "-ac -screen 0 1960x2000x24" .Rscript $@' >/root/bin/R_3.6.3/bin/Rscript && \
        chmod +x /root/bin/R_3.6.3/bin/Rscript

        # output
        OUTPUT_NAME="$(basename ~{querySketch} | sed 's/.msh$//')_to_$(basename ~{referenceSketch} | sed 's/.msh$//').tbl"

        # generate table
        mash dist -p ~{threadCount} -t ~{extraArguments} ~{querySketch} ~{referenceSketch} > $OUTPUT_NAME

        # generate plot
        ln -s $OUTPUT_NAME combined.tbl
        head -n 1 combined.tbl | awk '{for (i=2; i <=NF; i++) print $i}' |awk -F "/" '{print $NF}' > key
        Rscript $PLOT_R ${OUTPUT_NAME%.tbl}

        # cleanup
        rm combined.tbl

	>>>
	output {
		File table = glob("*.tbl")
		File plot = glob("*.png")
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}


task mashScreen {
    input {
        File readFile
        String extraArguments=""
        Int memSizeGB = 24
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_mash:latest"
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

        OUTPUT=$(basename ~{readFile} | sed 's/.msh$//').txt
        mash screen -p ~{threadCount} -w $MASH_REFSEQ ~{extraArguments} ~{readFile} > $OUTPUT

	>>>
	output {
		File screenOut = glob("*.txt")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
    }
}

