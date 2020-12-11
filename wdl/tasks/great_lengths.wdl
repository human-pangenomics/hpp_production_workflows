version 1.0

import "extract_reads.wdl" as extractReads_t
import "shard_reads.wdl" as shardReads_t
import "arithmetic.wdl" as arithmetic_t

workflow runGreatLengths {

    input {
        Array[File] reads
        File? referenceFasta
        String identifier="sample"
        Int histogramMinLength=0
        Int histogramMaxLength=0
        Int fileExtractionDiskSizeGB = 256
        String dockerImage = "tpesout/hpp_great_lengths:latest"
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

    # run greatLengths
    scatter (readFile in readsExtracted.extractedRead) {
        call greatLengths {
            input:
                readFile=readFile,
                histogramMinLength=histogramMinLength,
                histogramMaxLength=histogramMaxLength,
                diskSizeGB=maxReadSize.value + 10,
                dockerImage=dockerImage
        }
    }

    # consolidate
    call consolidateGreatLengths {
        input:
            greatLengthsTarballs=greatLengths.outputTarball,
            greatLengthsReports=greatLengths.reportForConsolidation,
            identifier=identifier,

    }

	output {
		File greatLengthsTarball = consolidateGreatLengths.outputTarball
		File greatLengthsReport = consolidateGreatLengths.outputReport
	}
}


task greatLengths {
    input {
        File readFile
        Int histogramMinLength = 0
        Int histogramMaxLength = 0
        Int memSizeGB = 2
        Int threadCount = 2
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_great_lengths:latest"
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

        # localize file for index
        ln -s ~{readFile}
        FILE=$(basename ~{readFile})

        # get output name
        OUTPUT=$(basename ~{readFile} | sed 's/.f[aq]\(st[aq]\)*$//')

        # hist parameters
        if [[ ~{histogramMinLength} -eq 0 && ~{histogramMaxLength} -eq 0 ]] ; then
            HIST_PARAM="--hist_auto_bounds"
        else
            HIST_PARAM="--hist_min ~{histogramMinLength} --hist_max ~{histogramMaxLength}"
        fi

        # sketch
        great_lengths.py -i $FILE -o $OUTPUT $HIST_PARAM

        # for consolidation
        cat $OUTPUT/report.tsv | sed "s/^/$OUTPUT\t/" >$OUTPUT.report.tsv

        # rename and tar
        ls $OUTPUT/ | xargs -n 1 -I{} mv $OUTPUT/{} $OUTPUT/$OUTPUT.{}
        tar czvf $OUTPUT.tar.gz $OUTPUT/

	>>>
	output {
		File outputTarball = glob("*.tar.gz")[0]
		File reportForConsolidation = glob("*.report.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task consolidateGreatLengths {
    input {
        Array[File] greatLengthsTarballs
        Array[File] greatLengthsReports
        String identifier="sample"
        Int memSizeGB = 4
        Int threadCount = 4
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_great_lengths:latest"
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

        # extract all tarballs
        mkdir ~{identifier}_great_lengths/
        cd ~{identifier}_great_lengths/
        for tb in ~{sep=" " greatLengthsTarballs} ; do
            tar xvf $tb &
        done
        wait
        cd ..
        tar czvf ~{identifier}_great_lengths.tar.gz ~{identifier}_great_lengths/

        # concat all reports
        for rp in ~{sep=" " greatLengthsReports} ; do
            cat $rp >>~{identifier}.report.tsv
        done
	>>>
	output {
		File outputTarball = glob("*.tar.gz")[0]
		File outputReport = glob("*.report.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

