version 1.0

workflow runBamCoverage{
    call bamCoverage
    output{
        File counts = bamCoverage.counts
        File coverageGz = bamCoverage.coverageGz
        Float coverageMeanFloat = bamCoverage.coverageMean
        Float coverageSdFloat = bamCoverage.coverageSD
    }
}

task bamCoverage{
    input{
        String sampleName
        String sampleSuffix
        String platform
        Array[File] bamFiles
        Int minMAPQ
        File assemblyFastaGz
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=256
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
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
        
        # Extract assembly and index
        FILENAME=`basename ~{assemblyFastaGz}`
        PREFIX="${FILENAME%.*.*}"
        gunzip -c ~{assemblyFastaGz} > $PREFIX.fa
        samtools faidx $PREFIX.fa

        if (( ~{length(bamFiles)} == 1 ))
        then 
            samtools depth -aa -Q ~{minMAPQ - 1} ~{sep=" " bamFiles}  > ~{sampleName}.depth
        else
            samtools depth -aa -Q ~{minMAPQ - 1} ~{sep=" " bamFiles} | awk '{sum=0; for (i=3; i<=NF; i++) { sum+= $i } {print $1,$2,sum}}' > ~{sampleName}.depth
        fi
        # Convert the output of samtools depth into a compressed format
        ${DEPTH2COV_C} -d ~{sampleName}.depth -f ${PREFIX}.fa.fai -o ~{sampleName}.~{sampleSuffix}.~{platform}.cov
        # Convert cov to counts
        ${COV2COUNTS_C} -i ~{sampleName}.~{sampleSuffix}.~{platform}.cov -o ~{sampleName}.~{sampleSuffix}.~{platform}.counts
        # Calculate mean and standard deviation
        python3 ${CALC_MEAN_SD_PY} --countsInput ~{sampleName}.~{sampleSuffix}.~{platform}.counts --meanOutput cov_mean.txt --sdOutput cov_sd.txt
        gzip ~{sampleName}.~{sampleSuffix}.~{platform}.cov
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File coverageGz = glob("*.cov.gz")[0]
        File counts = glob("*.counts")[0]
        Float coverageMean = read_float("cov_mean.txt") 
        Float coverageSD = read_float("cov_sd.txt")
    }
}

