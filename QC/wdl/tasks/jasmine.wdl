version 1.0

# This is a task level wdl workflow to run Jasmine

workflow runJasmine {
    input{
        File HiFiFile
        File OntFile
        String FileList
        String SV_like_errors
        String dockerImage

    }

    call Jasmine{
        input:
            HiFiFile = HiFiFile,
            OntFile = OntFile,
            FileList = FileList,
            SV_like_errors = SV_like_errors,
            dockerImage = dockerImage
    }
    output{
        File outputFile = Jasmine.outputFile
    }
}

task Jasmine{
    input{
        File HiFiFile
        File OntFile
        String FileList = "SV_filelist.txt"
        String SV_like_errors
        String dockerImage = "quay.io/biocontainers/jasminesv:1.1.4--hdfd78af_0"

    }
    # jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 --output_genotypes file_list="${OUTPUT_DIR}"/SV_filelist.txt out_file="${OUTPUT_DIR}"/"SV_like_errors.vcf"
    
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

        echo ~{OntFile} >> ~{FileList}
        echo ~{HiFiFile} >> ~{FileList}

        jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 --output_genotypes file_list=~{FileList} out_file=~{SV_like_errors}

    >>>
    output{
        File outputFile = glob("SV_like_errors.vcf")[0]
    }
    runtime{
        docker: dockerImage
    }
}
