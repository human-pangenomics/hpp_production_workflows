version 1.0

workflow runExtractReads {
    input {
        Array[File] inputFiles
        File? referenceFile
        Int threadCount=1
        String dockerImage
    }

    scatter (file in inputFiles) {
        call extractReads {
            input:
                readFile=file,
                referenceFile=referenceFile,
                threadCount=threadCount,
                dockerImage=dockerImage
        }
    }

    output {
        Array[File] reads = extractReads.extractedRead
    }
}

task extractReads {
    input {
        File readFile
        File? referenceFile
        Int threadCount=1
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
        # load samtools
        export PATH=$PATH:/root/bin/samtools_1.9/

        FILENAME=$(basename -- "~{readFile}")
        PREFIX="${FILENAME%.*}"
        SUFFIX="${FILENAME##*.}"

        mkdir output

        if [[ "$SUFFIX" == "bam" ]] ; then
            samtools fastq -@~{threadCount} ~{readFile} > output/${PREFIX}.fq
        elif [[ "$SUFFIX" == "cram" ]] ; then
            if [[ ! -f "~{referenceFile}" ]] ; then
                echo "Could not extract $FILENAME, reference file not supplied"
                exit 1
            fi
            ln -s ~{referenceFile}
            samtools fastq -@~{threadCount} --reference `basename ~{referenceFile}` ~{readFile} > output/${PREFIX}.fq
        elif [[ "$SUFFIX" == "gz" ]] ; then
            gunzip -k -c ~{readFile} > output/${PREFIX}
        elif [[ "$SUFFIX" != "fastq" ]] && [[ "$SUFFIX" != "fq" ]] && [[ "$SUFFIX" != "fasta" ]] && [[ "$SUFFIX" != "fa" ]] ; then
            echo "Unsupported file type: ${SUFFIX}"
            exit 1
        fi
    >>>

    output {
        File extractedRead = flatten([glob("output/*"), [readFile]])[0]
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
    }

    parameter_meta {
        readFile: {description: "Reads file in BAM, FASTQ, or FASTA format (optionally gzipped)"}
    }
}