version 1.0

workflow runFillGapsMale {
    input {
        File fai
        File autosome_nonCntr_bed
        File sex_nonCntr_bed
        File autosome_cntr_bed
        File sex_cntr_bed 
    }
    call fillGaps {
        input:
            fai = fai,
            autosome_nonCntr_bed = autosome_nonCntr_bed,
            sex_nonCntr_bed = sex_nonCntr_bed,
            autosome_cntr_bed = autosome_cntr_bed,
            sex_cntr_bed = sex_cntr_bed
    }
    output {
       File autosome_nonCntr_filled_bed = fillGaps.autosome_nonCntr_filled_bed
       File sex_nonCntr_filled_bed = fillGaps.sex_nonCntr_filled_bed
       File autosome_cntr_filled_bed = fillGaps.autosome_cntr_filled_bed
       File sex_cntr_filled_bed = fillGaps.sex_cntr_filled_bed
    }
}

task fillGaps {
    input {
        File fai
        File autosome_nonCntr_bed
        File sex_nonCntr_bed
        File autosome_cntr_bed
        File sex_cntr_bed
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=16
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
        
        # find and assign gaps (will produce bed files with the suffix "gaps.bed")
        python3 ${ASSIGN_GAPS_PY} --fai ~{fai} ~{autosome_nonCntr_bed} ~{sex_nonCntr_bed} ~{autosome_cntr_bed} ~{sex_cntr_bed}
       

        FILENAME_1=$(basename ~{autosome_nonCntr_bed})
        FILENAME_2=$(basename ~{sex_nonCntr_bed})
        FILENAME_3=$(basename ~{autosome_cntr_bed})
        FILENAME_4=$(basename ~{sex_cntr_bed})

        PREFIX_1=${FILENAME_1%.bed}
        PREFIX_2=${FILENAME_2%.bed}
        PREFIX_3=${FILENAME_3%.bed}
        PREFIX_4=${FILENAME_4%.bed}
        
        cat ~{autosome_nonCntr_bed} ${PREFIX_1}.gaps.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX_1}.filled_gaps.bed
        cat ~{sex_nonCntr_bed} ${PREFIX_2}.gaps.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX_2}.filled_gaps.bed
        cat ~{autosome_cntr_bed} ${PREFIX_3}.gaps.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX_3}.filled_gaps.bed
        cat ~{sex_cntr_bed} ${PREFIX_4}.gaps.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX_4}.filled_gaps.bed
           
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
    }
    output {
        File autosome_nonCntr_filled_bed = glob(basename("${autosome_nonCntr_bed}", ".bed") + ".filled_gaps.bed")[0]
        File sex_nonCntr_filled_bed = glob(basename("${sex_nonCntr_bed}", ".bed") + ".filled_gaps.bed")[0]
        File autosome_cntr_filled_bed = glob(basename("${autosome_cntr_bed}", ".bed") + ".filled_gaps.bed")[0]
        File sex_cntr_filled_bed = glob(basename("${sex_cntr_bed}", ".bed") + ".filled_gaps.bed")[0]
    }
}

