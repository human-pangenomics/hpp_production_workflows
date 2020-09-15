version 1.0

workflow runDipcall {
	call dipcall
}

task dipcall {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File referenceFasta
        Boolean isMaleSample
        Boolean referenceIsHS38 = true
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
        PATH="/root/bin/samtools_1.9:$PATH"

        # get output base
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/.[pm]at$//')
        mkdir dipcall_$PREFIX

        # prep fastas
        ln -s ~{referenceFasta}
        samtools faidx `basename ~{referenceFasta}`

        # initialize script
        cmd=( /opt/dipcall/dipcall.kit/run-dipcall )

        # male samples need PAR region excluded
        if [[ ~{isMaleSample} == true ]]; then
            if [[ ~{referenceIsHS38} ]]; then
                cmd+=( -x /opt/dipcall/dipcall.kit/hs38.PAR.bed )
            else
                cmd+=( -x /opt/dipcall/dipcall.kit/hs37d5.PAR.bed )
            fi

        fi

        # finalize script
        cmd+=( dipcall_$PREFIX/$PREFIX )
        cmd+=( `basename ~{referenceFasta}` )
        cmd+=( ~{assemblyFastaPat} )
        cmd+=( ~{assemblyFastaMat} )

        # generate makefile
        "${cmd[@]}" >$PREFIX.mak

        # run dipcall
        make -j 2 -f $PREFIX.mak

        # finalize
        tar czvf $PREFIX.dipcall.tar.gz dipcall_$PREFIX/
        cp dipcall_$PREFIX/$PREFIX.dip.bed $PREFIX.dipcall.bed
        cp dipcall_$PREFIX/$PREFIX.dip.vcf.gz $PREFIX.dipcall.vcf.gz

	>>>
	output {
		File outputTarball = glob("*.dipcall.tar.gz")[0]
		File outputVCF = glob("*.dipcall.vcf.gz")[0]
		File outputBED = glob("*.dipcall.bed")[0]
	}
    runtime {
        cpu: 2
        memory: "4 GB"
        docker: dockerImage
    }
}
