version 1.0 

workflow asset{
    input {
        String sampleName
        Array[File] ontPafFiles
        Array[File] hifiPafFiles
        Array[File] hicBamFiles
        File bionanoAlignmentTarGz
        File assembly
        Float ontCoverageMean
        Float ontCoverageSD
        Float hifiCoverageMean
        Float hifiCoverageSD
        Int preemptible=2
    }
    call ast_pbTask as ontAssetTask{
        input:
            sampleName = "${sampleName}.ont",
            pafFiles = ontPafFiles,
            coverageMean = ontCoverageMean,
            coverageSD = ontCoverageSD,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = preemptible
    }
    call ast_pbTask as hifiAssetTask{
        input:
            sampleName = "${sampleName}.hifi",
            pafFiles = hifiPafFiles,
            coverageMean = hifiCoverageMean,
            coverageSD = hifiCoverageSD,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = preemptible
    }
    call ast_hicTask as hicAssetTask{
        input:
            sampleName = sampleName,
            bamFiles = hicBamFiles,
            assembly = assembly,
            memSize = 32,
            threadCount = 8,
            diskSize = 512,
            preemptible = preemptible
    }
    call ast_bionTask as bionanoAssetTask{
        input:
            sampleName = sampleName,
            alignmentTarGz = bionanoAlignmentTarGz,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = preemptible
    }
    output{
        File hicBed = hicAssetTask.supportBed
        File hifiBed = hifiAssetTask.supportBed
        File ontBed = ontAssetTask.supportBed
        File bionanoBed = bionanoAssetTask.supportBed
        File gapsBed = hicAssetTask.gapsBed
        File hicCoverage = hicAssetTask.coverage
    }
}

task ast_bionTask{
    input{
        String sampleName
        File alignmentTarGz
        Int minCov=10
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
        Int preemptible
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

        tar -xzf ~{alignmentTarGz} --strip-components 1
        rmap=*_r.cmap
        qmap=*_q.cmap
        xmap=*.xmap
        keyfn=*_key.txt
        ast_bion_bnx -m~{minCov} $rmap $qmap $xmap $keyfn > ~{sampleName}.bionano.bed
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File supportBed = "~{sampleName}.bionano.bed"
    }

}



task ast_hicTask{
    input{
        String sampleName
        Array[File] bamFiles
        File assembly
        Int minMAPQ=20
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
        Int preemptible
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

        detgaps <(zcat ~{assembly}) > ~{sampleName}.gaps.bed
        ast_hic -q ~{minMAPQ} ~{sampleName}.gaps.bed ~{sep=" " bamFiles} > ast_hic.bed 2> ast_hic.log
        mv HC.base.cov ~{sampleName}.hic.cov
        
        GENOME_SIZE=`samtools view -H  ~{bamFiles[0]} | awk '{if($1 == "@SQ") {sum += substr($3,4,length($3))}} END {print sum}'`
        # Calculate mean coverage of aligned reads
        cat ~{sampleName}.hic.cov | awk -v GENOME_SIZE="$GENOME_SIZE" '{if( substr($1,1,1) != ">" ) {sum_cov += $3 * ($2 - $1 + 1)}} END {printf "%.2f\n", sqrt(sum_cov/GENOME_SIZE)}' > cov_mean.txt
        # Calculate standard deviation of base-level coverages (Considering coverages below 2.5 * COVERAGE_MEAN)
        cat ~{sampleName}.hic.cov | awk -v COVERAGE_MEAN=`cat cov_mean.txt` -v GENOME_SIZE="$GENOME_SIZE" '{if (( substr($1,1,1) != ">" ) && ($3 < 2.5 * COVERAGE_MEAN)) {sum += ($2 - $1 + 1) * (($3 - COVERAGE_MEAN)^2)}} END {printf "%.2f\n", sqrt(sum/GENOME_SIZE)}' > cov_sd.txt
        # calculate the min coverage threshold for asset
        # using the formula, min( max(10, mean - 2 x sd), 20)
        MIN_COVERAGE_ASSET=`awk -v mean=$(cat cov_mean.txt) -v sd=$(cat cov_sd.txt) 'BEGIN {min_cov = mean - 2 * sd; if (min_cov < 10) {min_cov=10}; if (min_cov > 20) {min_cov=20}; printf "%d",min_cov}'`
        cat ~{sampleName}.hic.cov | awk -v min_cov="${MIN_COVERAGE_ASSET}" '{if(substr($1,1,1) == ">"){chr=substr($1,2,5)} else if ($3 >= min_cov) {print chr"\t"$1"\t"$2}}' > ~{sampleName}.hic.bed
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File supportBed = "~{sampleName}.hic.bed"
        File gapsBed = "~{sampleName}.gaps.bed"
        File coverage = "~{sampleName}.hic.cov"
    }

}

task ast_pbTask{
    input{
        String sampleName
        Array[File] pafFiles
        Float coverageMean
        Float coverageSD
        # runtime configurations
        Int memSize
        Int threadCount
        Int diskSize
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
        Int preemptible
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
        
        # calculate the min coverage threshold for asset
        # using the formula, max(5, mean - 3 x sd)
        MIN_COVERAGE_ASSET=`awk -v mean=~{coverageMean} -v sd=~{coverageSD} 'BEGIN {min_cov = mean - 3 * sd; if (min_cov < 5) {min_cov=5}; if (min_cov > 10) {min_cov=10}; printf "%d",min_cov}'`
        # calculate the min coverage threshold for asset
        # using the formula, max(2.5 * mean, mean + 3 x sd)
        MAX_COVERAGE_ASSET=`awk -v mean=~{coverageMean} -v sd=~{coverageSD} 'BEGIN {max_cov = mean + 3 * sd; if (max_cov < (2.5 * mean)) {max_cov=2.5 * mean}; printf "%d",max_cov}'`
        # run asset to find supportive blocks
        ast_pb -m${MIN_COVERAGE_ASSET} -M${MAX_COVERAGE_ASSET} ~{sep=" " pafFiles} > ~{sampleName}.bed

    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output{
        File supportBed = "~{sampleName}.bed"
    }
}
