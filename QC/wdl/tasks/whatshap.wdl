version 1.0

workflow runWhatshapStats {

    input {
        Array[File] queryVcfs
        File truthVcf
        String outputIdentifier="sample"
    }

    scatter (queryVcf in queryVcfs) {
        call whatshapAnalysis {
            input:
                queryVcf=queryVcf,
                truthVcf=truthVcf
        }
    }

    call coalesceResults {
        input:
            tarballs=whatshapAnalysis.outputTarball,
            outputIdentifier=outputIdentifier
    }

    output {
        File outputTarball = coalesceResults.outputTarball
    }
}

task whatshapAnalysis {
    input {
        File queryVcf
        File truthVcf
        Int threadCount = 4
        Int memoryGB = 4
        String dockerImage="tpesout/whatshap:latest"
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

        OUTPUT_BASE=`basename ~{queryVcf} | sed 's/.gz$//' | sed 's/.vcf$//' | sed 's/\./_/g'`
        OUTPUT_DIR=${OUTPUT_BASE}.wh_stats
        OUTPUT=${OUTPUT_DIR}/${OUTPUT_BASE}
        mkdir $OUTPUT_DIR

        QUERY=`basename ~{queryVcf} | sed 's/.gz$//'`
        TRUTH=`basename ~{truthVcf} | sed 's/.gz$//'`
        if [[ "~{queryVcf}" =~ ".gz" ]] ; then
            zcat ~{queryVcf} | sed 's|^\(#[^#].*\t\).*|\1sample|' >$QUERY
        else
            cat ~{queryVcf} | sed 's|^\(#[^#].*\t\).*|\1sample|' >$QUERY
        fi
        if [[ "~{truthVcf}" =~ ".gz" ]] ; then
            zcat ~{truthVcf} | sed 's|^\(#[^#].*\t\).*|\1sample|' >$TRUTH
        else
            cat ~{truthVcf} | sed 's|^\(#[^#].*\t\).*|\1sample|' >$TRUTH
        fi


        # get stats
        whatshap stats \
            --tsv ${OUTPUT}.stats.tsv \
            --block-list ${OUTPUT}.blocks.txt \
            --chr-lengths /root/chr_lengths \
            $QUERY \
            >${OUTPUT}.stats.txt

        # get comparison
        whatshap compare \
            --tsv-pairwise ${OUTPUT}.pairwise.tsv \
            --tsv-multiway ${OUTPUT}.multiway.tsv \
            --switch-error-bed ${OUTPUT}.switch_error.bed \
            --plot-blocksizes ${OUTPUT}.blocksizes.pdf \
            --plot-sum-of-blocksizes ${OUTPUT}.sum_of_blocksizes.pdf \
            --longest-block-tsv ${OUTPUT}.longest_block.tsv \
            $QUERY $TRUTH \
            >$OUTPUT.compare.txt

        # tarball it
        tar czvf ${OUTPUT_DIR}.tar.gz ${OUTPUT_DIR}/


    >>>

    output {
        File outputTarball = glob("*.tar.gz")[0]
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memoryGB + "GB"
    }
}

task coalesceResults {
    input {
        Array[File] tarballs
        String outputIdentifier
        Int threadCount = 4
        Int memoryGB = 4
        String dockerImage="tpesout/whatshap:latest"
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

        for TB in ~{sep=" " tarballs} ; do
            tar xvf $TB &
        done
        wait

        # this ugly sed basically basenames the vcfs with termination ".vcf" or ".vcf.gz"
        for DIR in `ls | grep wh_stats` ; do
            if [[ ! -f "~{outputIdentifier}.stats.tsv" ]] ; then
                cat $DIR/*stats.tsv | sed 's|\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)\t|\t\1\t|' > ~{outputIdentifier}.stats.tsv
                cat $DIR/*pairwise.tsv | sed 's|\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)|\t\1\t\3|' > ~{outputIdentifier}.pairwise.tsv
            else
                cat $DIR/*stats.tsv | sed 's|\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)\t|\t\1\t|' | grep -v "^#" >> ~{outputIdentifier}.stats.tsv
                cat $DIR/*pairwise.tsv | sed 's|\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)\t/.*/\(.*\.vcf\(.gz\)\{0,1\}\)|\t\1\t\3|' | grep -v "^#" >> ~{outputIdentifier}.pairwise.tsv
            fi
        done

        # merge into one consistent file
        python3 /root/merge_results.py ~{outputIdentifier}.stats.tsv ~{outputIdentifier}.pairwise.tsv >~{outputIdentifier}.full.tsv

        # tarball
        tar czvf ~{outputIdentifier}.wh_stats.tar.gz *.wh_stats ~{outputIdentifier}.stats.tsv ~{outputIdentifier}.pairwise.tsv ~{outputIdentifier}.full.tsv
    >>>

    output {
        File outputTarball = outputIdentifier + ".wh_stats.tar.gz"
        File statsOutput = outputIdentifier + ".stats.tsv"
        File pairwiseOutput = outputIdentifier + ".pairwise.tsv"
        File fullOutput = outputIdentifier + ".full.tsv"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memoryGB + "GB"
    }
}
