version 1.0

workflow runPepperMarginDeepVariant{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        ## Model types can be:
        ##     "ont_r9_guppy5_sup"
        ##     "ont_r10_q20"
        ##     "hifi"
        String pmdvModelType
        Int minMAPQ = 0
        String includeSupplementary="True"
    }
    call pmdv{
        input:
            modelType = pmdvModelType,
            assemblyFastaGz = assemblyFastaGz,
            bam = bam,
            bamIndex = bamIndex,
            includeSupplementary = includeSupplementary,
            minMAPQ = minMAPQ,
            diskSize= 4 * ceil(size(bam, "GB")) + 64
    }
    output{
        File vcfGz = pmdv.vcfGz
    }
}

task pmdv{
    input{
        File bam
        File bamIndex
        File assemblyFastaGz
        Int minMAPQ
        String includeSupplementary="False"
        String modelType
        # runtime configurations
        Int memSize=256
        Int threadCount=64
        Int diskSize=512
        String dockerImage="kishwars/pepper_deepvariant:r0.6"
        Int preemptible=2
        String zones="us-west2-a"
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
        
        ## hard link the bam file to the working directory and produce its index file
        BAM_NAME=$(basename ~{bam})
        BAM_PREFIX=${BAM_NAME%%.bam}
        ln -f ~{bam} > ${BAM_PREFIX}.bam
        ln -f ~{bamIndex} > ${BAM_PREFIX}.bam.bai
        # Update the creation time of the index file
        touch -c ${BAM_PREFIX}.bam.bai 

        ## unzip the fasta file and produce its index
        gunzip -c ~{assemblyFastaGz} > asm.fa
        samtools faidx asm.fa

        MORE_OPTIONS="--dv_min_mapping_quality ~{minMAPQ} --pepper_min_mapping_quality ~{minMAPQ}"
        if [ ~{includeSupplementary} == "True" ]; then
            MORE_OPTIONS="${MORE_OPTIONS} --pepper_include_supplementary"
        fi

        mkdir output

        ## call pepper-margin-deepvariant
        run_pepper_margin_deepvariant call_variant \
        -b  ${BAM_PREFIX}.bam \
        -f  asm.fa \
        -o  output \
        -t ~{threadCount} \
        --~{modelType} ${MORE_OPTIONS}

        mv output/PEPPER_MARGIN_DEEPVARIANT_OUTPUT.vcf.gz output/${BAM_PREFIX}.vcf.gz
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File vcfGz = glob("output/*.vcf.gz")[0]
    }
}

