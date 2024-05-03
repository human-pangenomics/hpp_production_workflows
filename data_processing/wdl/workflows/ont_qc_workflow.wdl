version 1.0

import "../tasks/calc_ont_summary_stats.wdl" as calc_ont_stats_wf
import "../tasks/ntsm.wdl" as ntsm_check

workflow ont_qc_wf {
    input {
        File ont_reads
        File ont_sequencing_summary_files
        Array[File] other_reads
        String file_id
    }

    parameter_meta {
        ont_reads: "Reads to QC. Should be single file. Can be .fq.gz, fastq.gz, or bam."
        ont_sequencing_summary_files: "Sequencing summary file output by Dorado/Guppy."
        other_reads: "Trusted reads to compare against to check for sample swaps. Can be fastq, fastq.gz, or bam."
        file_id: "Id to use for output naming."
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
    }
    
    ## aggregate coverage information (broken down by read length)
    call calc_ont_stats_wf.run_calc_ont_stats as ont_stats {
        input: 
            sequencing_summary = ont_sequencing_summary_files
    }

    ## check for sample swaps
    call ntsm_check.ntsm_workflow as ntsm_wf {
        input:
            input_reads_1  = [ont_reads],
            input_reads_2  = other_reads,
            sample_id      = file_id,
            read_1_type    = "ont",
            read_2_type    = "other"
    }

    ## aggregate results
    call summarize_ont_qc {
        input:
            ntsm_output         = ntsm_wf.ntsm_eval_out,
            pass_summary_stats  = ont_stats.pass_summary_stats ,
            fail_summary_stats  = ont_stats.fail_summary_stats,
            file_name           = file_id
    }


    output {
        File ont_qc_summary       = summarize_ont_qc.summary_file

        ## ont summary stats
        File ont_pass_stats       = ont_stats.pass_summary_stats 
        File ont_fail_stats       = ont_stats.fail_summary_stats

        ## ntsm output (check for sample swaps)
        File ont_ntsm_counts      = ntsm_wf.ntsm_count_1
        File ext_ntsm_counts      = ntsm_wf.ntsm_count_2
        File ntsm_eval            = ntsm_wf.ntsm_eval_out
    }
}


task summarize_ont_qc {
    input {
        File ntsm_output
        File pass_summary_stats
        File fail_summary_stats
        String file_name

        Int memSizeGB   = 4
        Int threadCount = 1
        Int disk_dize   = 16        
        Int preempts    = 2
    }

    command <<<
        set -eux -o pipefail

        # Define file paths
        pass_summary="~{pass_summary_stats}"
        fail_summary="~{fail_summary_stats}"
        nstm_file="~{ntsm_output}"
        filename="~{file_name}"

        # Create output file with a header
        echo -e "file_name\tpass_read_N50\tpass_Gb\tpass_coverage\tpass_100kb+\tpass_200kb+\tpass_300kb+\tpass_400kb+\tpass_500kb+\tpass_1Mb+\tpass_whales\tfail_read_N50\tfail_Gb\tfail_coverage\tfail_100kb+\tfail_200kb+\tfail_300kb+\tfail_400kb+\tfail_500kb+\tfail_1Mb+\tfail_whales\tntsm_score\tntsm_result" > aggregated_output.tsv

        # Extract data from pass_summary including the file name and skip header
        # Note that the file name is formatted as ['file_name_pass.txt']
        # we want file_name
        pass_data=$(awk -F'\t' 'NR > 1 {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' $pass_summary)

        # Extract data from fail_summary, skip the file name and header
        fail_data=$(awk -F'\t' 'NR > 1 {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' $fail_summary)

        # Extract relevant columns from nstm_file and skip any header if present
        nstm_data=$(awk -F'\t' '{print $3"\t"$4}' $nstm_file)

        # Append data to the output file
        echo -e "$filename\t$pass_data\t$fail_data\t$nstm_data" >> "~{file_name}.summary.tsv"
    >>>

    output {
        File summary_file = "~{file_name}.summary.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + disk_dize + " SSD"
        docker: "humanpangenomics/ont_summary_stats@sha256:ec6eba759f8f91cf461262ba5155aa951b77acf241f764dda4b58f7d5247d96b"
        preemptible: preempts
    }
}


