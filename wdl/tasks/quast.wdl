version 1.0

workflow runQuast {
	call quast
}

task quast {
    input {
        File assemblyFasta
        File? referenceFasta
        String extraArguments="--min-identity 80 --fragmented --large"
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

        # initilization
        ln -s ~{assemblyFasta}
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $file =~ \.gz$ ]]; then
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # init quast command
        cmd=(python /root/tools/quast/quast-5.0.2/quast-lg.py )
        cmd+=( -t ~{threadCount} )
        cmd+=( -o $PREFIX.quast )

        # include reference fasta if supplied
        if [[ -f "~{referenceFasta}" ]]; then
            ln -s ~{referenceFasta}
            REF_FILENAME=$(basename -- "~{assemblyFasta}")
            if [[ $file =~ \.gz$ ]]; then
                gunzip $REF_FILENAME
                REF_FILENAME="${REF_FILENAME%.gz}"
            fi
            cmd+=( -r $REF_FILENAME )
        fi

        # include extra arguments if supplied
        if [[ ! -z "~{extraArguments}" ]] ; then
            cmd+=( ~{extraArguments} )
        fi

        # finalize command
        cmd+=( $ASM_FILENAME )

        # run command
        "${cmd[@]}"

        # save output
        tar czvf $PREFIX.quast.tar.gz $PREFIX.quast

	>>>
	output {
		File outputTarball = glob("*.quast.tar.gz")[0]
	}
    runtime {
        cpu: threadCount
        memory: "8 GB"
        docker: dockerImage
    }
}
