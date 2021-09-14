version 1.0 

workflow runCov2CountsContigWise{
    call cov2countsContigWise
    output {
        File contigCountsTarGz = cov2countsContigWise.contigCountsTarGz
        File contigCovsTarGz = cov2countsContigWise.contigCovsTarGz
    }
}

task cov2countsContigWise {
    input {
        File coverageGz
        # runtime configurations
        Int memSize=8
        Int threadCount=8
        Int diskSize=64
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
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

        FILENAME=$(basename ~{coverageGz})
        PREFIX=${FILENAME%.cov.gz}
        
        gunzip -c ~{coverageGz} > ${PREFIX}.cov
        mkdir covs counts
        # Make a separate cov file for each contig
        ${SPLIT_CONTIGS_COV_BIN} -c ${PREFIX}.cov -p covs/${PREFIX}
        # Count each contig-specific cov file
        for c in $(ls covs);do ${COV2COUNTS_BIN} -i covs/$c -o counts/${c/.cov/.counts}; echo $c" finished";done

        # Compress Counts files
        tar -cf ${PREFIX}.counts.tar counts
        gzip ${PREFIX}.counts.tar

        # Compress Cov files
        tar -cf ${PREFIX}.covs.tar covs
        gzip ${PREFIX}.covs.tar
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File contigCountsTarGz = glob("*.counts.tar.gz")[0]
        File contigCovsTarGz = glob("*.covs.tar.gz")[0]
    }
}

