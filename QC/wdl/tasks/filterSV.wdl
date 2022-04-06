version 1.0

# This is a task level wdl workflow to run Kishwar's python script to filter SV calls

workflow runFilterSV {
    input {
        File inputVcf
        String outputName
        String dockerImage

    }
    call Filter {
        input:
            inputVcf = inputVcf,
            outputName = outputName,
            dockerImage = dockerImage
    }
    output {
        File outputFile = Filter.outputFile
    }
}

task Filter{
    input {
        File inputVcf
        String outputName
        String dockerImage = "kishwars/t2t_polishing:0.1"
        # String dockerImage = "https://github.com/kishwarshafin/T2T_polishing_scripts/blob/master/docker/Dockerfile"

    }
    # python3 filter.py "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES}" > "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES_FILTERED}"
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
        ln -s /opt/filter.py .
        python3 filter.py ~{inputVcf} > ~{outputName}

    >>>
    output {
        File outputFile = glob("*filtered.vcf")[0]
    }
    runtime {
        docker: dockerImage
    }
}
