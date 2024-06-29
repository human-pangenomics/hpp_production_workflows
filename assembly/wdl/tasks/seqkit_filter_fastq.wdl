version 1.0

workflow seqkit_filter_fastq {

    call filter_fastq

    output {
        File seqkit_filtered_fq = filter_fastq.filteredFastq
    }
}


task filter_fastq {

    input {
        File input_fastq

        Int min_size    = -1
        Float min_q     = 10
        
        Int threadCount = 4
        Int memSizeGB   = 8
        Int addldisk    = 32

        String dockerImage = "quay.io/biocontainers/seqkit:0.15.0--0"
    }

    Int fastq_size = ceil(size(input_fastq, "GB"))    
    Int final_disk_dize = 2*fastq_size + addldisk

    String output_prefix = basename(input_fastq, ".gz")
    String output_name   = "~{output_prefix}.filtered.gz"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        seqkit seq \
            --min-len ~{min_size} \
            --min-qual ~{min_q} \
            ~{input_fastq} \
            -o ~{output_name}

    >>>

    output {

        File filteredFastq = output_name
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}