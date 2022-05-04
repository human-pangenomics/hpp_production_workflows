version 1.0

# This is a task level wdl workflow to run Jasmine

workflow runJasmine {
    input{
        Array[File] InputVCFs
        String SV_like_errors
        String? dockerImage
        Int? maxDist = 500
        Float? minSeqID = 0.3
        Int? specReads = 3
    }

    call Jasmine{
        input:
            InputVCFs = InputVCFs,
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
        Array[File] InputVCFs
        String SV_like_errors
        String dockerImage = "quay.io/biocontainers/jasminesv:1.1.4--hdfd78af_0"
        Int? maxDist = 500
        Float? minSeqID = 0.3
        Int? specReads = 3
        Int memSizeGB = 64
        Int threadCount = 32
        Int diskSizeGB = 64

    }
    
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        
        VCFS=(~{sep=' ' InputVCFs})

        for VCF in "${VCFS[@]}"
        do
                echo "$VCF" > SV_filelist.txt
        done

        jasmine max_dist=~{maxDist} min_seq_id=~{minSeqID} spec_reads=~{specReads} --output_genotypes file_list=SV_filelist.txt out_file=~{SV_like_errors}

    >>>
    output{
        File SV_filelist = glob("*.txt")[0]
        File outputFile = glob("SV_like_errors.vcf")[0]
    }
    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
