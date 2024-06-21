version 1.0


workflow concatVCF {

    call bcftoolsConcat

    output {
        File concatVcf=bcftoolsConcat.vcfOut
        File concatVcfIdx=bcftoolsConcat.vcfOutIdx
    }
}

task bcftoolsConcat {
    input {
        File vcf1
        File vcf2

        Int memSizeGB = 8
        Int threadCount = 4
        Int diskSizeGB = 50
        String dockerImage = "kishwars/pepper_deepvariant:r0.8"

    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        FILENAME1=$(basename -- "~{vcf1}")
        SUFFIX1="${FILENAME1##*.}"

        FILENAME2=$(basename -- "~{vcf2}")
        SUFFIX2="${FILENAME2##*.}"

        if [[ "$SUFFIX1" != "gz" ]] ; then
            bgzip -c ~{vcf1} > ./vcf1.vcf.gz
        else
            cp ~{vcf1} ./vcf1.vcf.gz
        fi

        if [[ "$SUFFIX2" != "gz" ]] ; then
            bgzip -c ~{vcf2} > ./vcf2.vcf.gz
        else
            cp ~{vcf2} ./vcf2.vcf.gz 
        fi

        VCF1_PREFIX=$(basename ~{vcf1})
        VCF2_PREFIX=$(basename ~{vcf2})

        tabix -p vcf ./vcf1.vcf.gz
        tabix -p vcf ./vcf2.vcf.gz

        mkdir output

        bcftools concat -a ./vcf1.vcf.gz ./vcf2.vcf.gz > output/${VCF1_PREFIX}_${VCF2_PREFIX}.vcf

        bgzip output/${VCF1_PREFIX}_${VCF2_PREFIX}.vcf
        tabix -p vcf output/${VCF1_PREFIX}_${VCF2_PREFIX}.vcf.gz
    >>>

    output {
        File vcfOut = glob("output/*.vcf.gz")[0]
        File vcfOutIdx = glob("output/*.vcf.gz.tbi")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 2
    }
}
