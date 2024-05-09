version 1.0

import "../tasks/fastqReadCoverage.wdl" as sum_fastq_bases
import "../tasks/ntsm.wdl" as ntsm_check

workflow hic_qc_wf {
    input {
        File hic_reads
        Array[File] other_reads
        String file_id
    }

    parameter_meta {
        hic_reads: "Reads to QC. Should be single file. Can be *.fq.gz or *.fastq.gz"
        other_reads: "Trusted reads to compare against to check for sample swaps. Can be fastq, fastq.gz, or bam."
        file_id: "Id to use for output naming."
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
    }
    
    ## sum bases in all reads (does not require that all reads are same length)
    call sum_fastq_bases.sumFastqReads as sum_reads {
        input: 
            inputFastq = [hic_reads],
            out_prefix = "~{file_id}_base_count.txt"
    }

    ## check for sample swaps
    call ntsm_check.ntsm_workflow as ntsm_wf {
        input:
            input_reads_1  = [hic_reads],
            input_reads_2  = other_reads,
            sample_id      = file_id,
            read_1_type    = "hic",
            read_2_type    = "other"
    }

    ## aggregate results
    call summarize_hic_qc {
        input:
            ntsm_output         = ntsm_wf.ntsm_eval_out,
            coverage_file       = sum_reads.coverageFile,
            file_name           = file_id
    }


    output {
        File hic_qc_summary     = summarize_hic_qc.summary_file

        ## ntsm output (check for sample swaps)
        File ont_ntsm_counts    = ntsm_wf.ntsm_count_1
        File ext_ntsm_counts    = ntsm_wf.ntsm_count_2
        File ntsm_eval          = ntsm_wf.ntsm_eval_out
    }
}


task summarize_hic_qc {
    input {
        File ntsm_output
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
        nstm_file="~{ntsm_output}"
        filename="~{file_name}"

        # Create output file with a header
        echo -e "file_name\ttotal_bp\tcoverage\tntsm_score\tntsm_result" > "~{file_name}.summary.tsv"


        # we want file_name
    
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

        # Extract relevant columns from nstm_file and skip any header if present
        nstm_data=$(awk -F'\t' '{print $3"\t"$4}' $nstm_file)
        
        # Append data to the output file
        echo -e "$filename\t$cov_line\t$nstm_data" >> "~{file_name}.summary.tsv"

        
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


