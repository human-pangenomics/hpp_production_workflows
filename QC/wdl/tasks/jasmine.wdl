version 1.0

# This is a task level wdl workflow to run Jasmine

workflow runJasmine {
    input{
        File HiFiFile
        File OntFile
        String SV_like_errors
        String? dockerImage
        Int? maxDist = 500
        Float? minSeqID = 0.3
        Int? specReads = 3
    }

    call Jasmine{
        input:
            HiFiFile = HiFiFile,
            OntFile = OntFile,
            SV_like_errors = SV_like_errors,
            dockerImage = dockerImage,
            maxDist = maxDist,
            minSeqID = minSeqID,
            specReads = specReads
    }
    output{
        File SV_filelist = Jasmine.SV_filelist
        File outputFile = Jasmine.outputFile
    }
}

task Jasmine{
    input{
        File HiFiFile
        File OntFile
        String SV_like_errors
        String dockerImage = "quay.io/biocontainers/jasminesv:1.1.4--hdfd78af_0"
        Int? maxDist = 500
        Float? minSeqID = 0.3
        Int? specReads = 3

    }
    # jasmine max_dist=500 min_seq_id=0.3 spec_reads=3 --output_genotypes file_list="${OUTPUT_DIR}"/SV_filelist.txt out_file="${OUTPUT_DIR}"/"SV_like_errors.vcf"
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        echo ~{OntFile} >> SV_filelist.txt
        echo ~{HiFiFile} >> SV_filelist.txt

        jasmine max_dist=~{maxDist} min_seq_id=~{minSeqID} spec_reads=~{specReads} --output_genotypes file_list=SV_filelist.txt out_file=~{SV_like_errors}

    >>>
    output{
        File SV_filelist = glob("*.txt")[0]
        File outputFile = glob("SV_like_errors.vcf")[0]
    }
    runtime{
        docker: dockerImage
    }
}
