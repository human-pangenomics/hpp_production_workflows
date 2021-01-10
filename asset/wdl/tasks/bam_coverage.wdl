version 1.0

workflow runBamCoverage{
    call bamCoverage
}

task bamCoverage{
    input{
        Array[File] bamFiles
        Int minMAPQ
        # runtime configurations
        Int memSize=16
        Int threadCount=4
        Int diskSize=256
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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

        if (( ~{length(bamFiles)} == 1 ))
        then 
            samtools depth -Q ~{minMAPQ - 1} ~{sep=" " bamFiles}  > coverage.bed
        else
            samtools depth -Q ~{minMAPQ - 1} ~{sep=" " bamFiles} | awk '{sum=0; for (i=3; i<=NF; i++) { sum+= $i } {print $1,$2,sum}}' > coverage.bed
        fi
        GENOME_SIZE=`samtools view -H  ~{bamFiles[0]} | awk '{if($1 == "@SQ") {sum += substr($3,4,length($3))}} END {print sum}'`
        # Calculate mean coverage of aligned reads
        cat coverage.bed | awk '{sum += $3} END {printf "%.2f\n", sum/NR}' > cov_mean.txt
        # Calculate standard deviation of base-level coverages (Considering coverages below 2.5 * COVERAGE_MEAN)
        cat coverage.bed | awk -v COVERAGE_MEAN=`cat cov_mean.txt` '{if ($3 < 2.5 * COVERAGE_MEAN) {sum += ($3 - COVERAGE_MEAN)^2}} END {printf "%.2f\n", sqrt(sum/NR)}' > cov_sd.txt
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        Float coverageMean = read_float("cov_mean.txt") 
        Float coverageSD = read_float("cov_sd.txt")
    }
}

