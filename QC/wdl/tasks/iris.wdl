version 1.0

# This is a task level wdl workflow to run IRIS to filter HIFI SV calls

workflow runIris {
    input{
        File genomeIn
        File readsIn
        File vcfIn
        String vcfOut
        String IrisOut
        String dockerImage

    }
    call Iris{
        input:
            genomeIn = genomeIn,
            readsIn = readsIn,
            vcfIn = vcfIn,
            vcfOut = vcfOut,
            IrisOut = IrisOut,
            dockerImage = dockerImage
    }
    output{
        File outputFile = Iris.outputFile
    }
}

task Iris{
    input{
        File genomeIn
        File readsIn
        File vcfIn
        String vcfOut
        String IrisOut
        String dockerImage = "quay.io/biocontainers/irissv:1.0.4--hdfd78af_2"

    }
    # iris --hifi --keep_long_variants --keep_files genome_in="${INPUT_DIR}"/chm13.draft_v0.9.chr20.fasta 
    # reads_in="${INPUT_DIR}"/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam vcf_in="${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED}" 
    # vcf_out="${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED_IRIS}"
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

        iris --hifi --keep_long_variants --keep_files genome_in=~{genomeIn} reads_in=~{readsIn} vcf_in=~{vcfIn} vcf_out=~{vcfOut} out_dir=~{IrisOut}

    >>>
    output{
        File outputFile = glob("*.iris.vcf")[0]
    }
    runtime{
        docker: dockerImage
    }
}
