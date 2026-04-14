version 1.0

import "../tasks/fastqReadCoverage.wdl" as sum_fastq_bases
import "../tasks/ntsm.wdl" as ntsm_check

workflow hic_qc_wf {
    input {
        Array[File] hic_reads
        String file_id
    }

    parameter_meta {
        hic_reads: "Paired-end reads to QC (R1 and R2). Can be *.fq.gz or *.fastq.gz"
        file_id: "Id to use for output naming."
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
    }
    
    ## sum bases in all reads (does not require that all reads are same length)
    call sum_fastq_bases.sumFastqReads as sum_reads {
        input:
            inputFastq = hic_reads,
            out_prefix = "~{file_id}_base_count.txt"
    }

    ## check for sample swaps
    call ntsm_check.ntsm_count as ntsm_wf {
        input:
            input_reads = hic_reads,
            sample_id   = file_id,
            read_type   = "hic"
    }

    ## aggregate results
    call summarize_hic_qc {
        input:
            coverage_file = sum_reads.coverageFile,
            file_name     = file_id
    }


    output {
        File hic_qc_summary  = summarize_hic_qc.summary_file

        ## ntsm count output (eval is run separately across all samples)
        File hic_ntsm_counts = ntsm_wf.ntsm_counts
    }
}


task summarize_hic_qc {
    input {
        File coverage_file
        String file_name

        Int memSizeGB   = 4
        Int threadCount = 1
        Int disk_dize   = 16
        Int preempts    = 2
    }

    command <<<
        set -eux -o pipefail

        # Define file paths
        coverage_txt="~{coverage_file}"
        filename="~{file_name}"

        # Create output file with a header
        echo -e "file_name\ttotal_bp\tcoverage" > "~{file_name}.summary.tsv"

        cov_line=$(awk '
            BEGIN {
              genome_size = 3.1e9;
              OFS="\t"
            }
            {
              total_bases = $1;
              coverage = (total_bases / genome_size);
              printf("%d\t%.2f\n", total_bases, coverage);
            }
            ' "$coverage_txt")

        # Append data to the output file
        echo -e "$filename\t$cov_line" >> "~{file_name}.summary.tsv"
    >>>

    output {
        File summary_file = "~{file_name}.summary.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + disk_dize + " SSD"
        docker: "registry.access.redhat.com/ubi7/ubi:latest"
        preemptible: preempts
    }
}


