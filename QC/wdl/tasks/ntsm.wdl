version 1.0

workflow ntsm_workflow {
    input {
        Array[File] input_fastqs_1
        Array[File] input_fastqs_2
        String sample_id   = "sample"
        String read_1_type = "type1"
        String read_2_type = "type2"
    }


    call ntsm_count as ntsm_count_1 {
        input:
            input_fastqs = input_fastqs_1,
            sample_id    = sample_id,
            read_type    = read_1_type
    }

    call ntsm_count as ntsm_count_2 {
        input:
            input_fastqs = input_fastqs_2,
            sample_id    = sample_id,
            read_type    = read_2_type
    }

    call ntsm_eval {
        input:
            count_1_file = ntsm_count_1.ntsm_counts,
            count_2_file = ntsm_count_2.ntsm_counts,
            sample_id    = sample_id,
            read_1_type  = read_1_type,
            read_2_type  = read_2_type
    }

   output {
       File ntsm_eval_out = ntsm_eval.ntsm_eval_out
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calls [ntsm](https://github.com/JustinChu/ntsm) to identify mismatched read files through k-mer analysis"
    }
}


task ntsm_count {
    input {
        Array[File] input_fastqs
        String sample_id
        String read_type

        Int count_kmer_size = 19

        Int memSizeGB   = 4
        Int threadCount = 4
        Int addldisk    = 10
        Int preempts    = 2
    }
    
    parameter_meta {
        input_fastqs: "Files must be in fastq.gz format"
        count_kmer_size: "k-mer size to use (sliding window is applied: highest count is stored)"
    }

    # Estimate disk size required
    Int input_fastq_size  = ceil(size(input_fastqs, "GB"))       
    Int final_disk_dize   = input_fastq_size * 2 + addldisk

    # Create output file name
    String output_counts_fn = "${sample_id}_${read_type}_counts.txt"

    command <<<
        set -eux -o pipefail

        ## localize kmers
        ln -s /opt/ntsm/data/31_AT.fa .
        ln -s /opt/ntsm/data/31_CG.fa .

        FASTQS=(~{sep=" " input_fastqs})

        ## Count number of k-mers found 
        /opt/ntsm/src/ntsmCount \
            -r 31_AT.fa \
            -a 31_CG.fa \
            -k ~{count_kmer_size} \
            -t ~{threadCount} \
            <(pigz -cd ${FASTQS}) \
            > ~{output_counts_fn}
    >>>

    output {
        File ntsm_counts = output_counts_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "humanpangenomics/ntsm:latest"
        preemptible: preempts
    }
}

task ntsm_eval {
    input {
        File count_1_file
        File count_2_file       
        String sample_id
        String read_1_type
        String read_2_type

        Int memSizeGB   = 4
        Int threadCount = 4
        Int diskSize    = 32
        Int preempts    = 2
    }

    # Create output file name
    String output_eval_fn = "${sample_id}_${read_1_type}_vs_${read_2_type}.txt"

    command <<<
        set -eux -o pipefail

        /opt/ntsm/src/ntsmEval \
            ~{count_1_file} \
            ~{count_2_file} \
            > ~{output_eval_fn}
    >>>

    output {
        File ntsm_eval_out = output_eval_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "humanpangenomics/ntsm:latest"
        preemptible: preempts
    }
}