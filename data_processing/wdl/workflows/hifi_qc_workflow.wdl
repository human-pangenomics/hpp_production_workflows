version 1.0

import "../tasks/methylationCheck.wdl" as meth_check
import "../tasks/read_stats.wdl" as read_stats_wf
import "../tasks/ntsm.wdl" as ntsm_check

workflow hifi_qc_wf {
    input {
        File hifi_reads
        Array[File] other_reads
        String sample_id
        Boolean perform_methylation_check = true  # Optional input to perform methylation check
    }

    parameter_meta {
        hifi_reads: "Reads to QC. Can be fastq, fastq.gz, or bam."
        other_reads: "Trusted reads to compare against to check for sample swaps. Can be fastq, fastq.gz, or bam."
        sample_id: "Id to use for output naming."
        perform_methylation_check: "Set to true to perform methylation check on BAM files."
    }

    meta {
        maintainer: "Andrew Blair"
        email: "apblair@ucsc.edu"
    }

    # Optional methylation check
    if (perform_methylation_check) {
        call meth_check.methylationWorkflow as methylation {
            input:
                bam_files = [hifi_reads],
                sample_id = sample_id
        }
    }

    call read_stats_wf.runReadStats as read_stats {
        input: 
            reads = [hifi_reads],
            identifier = sample_id
    }

    call ntsm_check.ntsm_workflow as ntsm_wf {
        input:
            input_reads_1 = [hifi_reads],
            input_reads_2 = other_reads,
            sample_id = sample_id,
            read_1_type = "hifi",
            read_2_type = "other"
    }

    call summarize_hifi_qc {
        input:
            ntsm_output = ntsm_wf.ntsm_eval_out,
            readstat_report = read_stats.ReadStatsReport,
            file_name = sample_id,
            methylation_report = select_first([if (perform_methylation_check) then methylation.methylation_summary else "/dev/null"])
    }

    output {
        ## readstats output
        File hifi_readstat_tar = read_stats.ReadStatsTarball 
        File hifi_readstat_report = read_stats.ReadStatsReport
        
        ## ntsm output
        File hifi_ntsm_counts = ntsm_wf.ntsm_count_1
        File ext_ntsm_counts = ntsm_wf.ntsm_count_2
        File ntsm_eval = ntsm_wf.ntsm_eval_out
        ## methylation check output
        File? methylation_report = select_first([if (perform_methylation_check) then methylation.methylation_summary else "/dev/null"])
        File hifi_qc_summary = summarize_hifi_qc.summary_file

  
    }
}

task summarize_hifi_qc {
    input {
        File ntsm_output
        File readstat_report
        String file_name
        File? methylation_report = "/dev/null"

        Int memSizeGB = 4
        Int threadCount = 1
        Int disk_size = 16        
        Int preempts = 2
    }

    command <<<
        set -eux -o pipefail

        # Log the input files
        echo "Read stats file: ~{readstat_report}"
        echo "NTSM output file: ~{ntsm_output}"
        echo "Methylation report: ~{methylation_report}"

        # Define file paths
        readstat_file="~{readstat_report}"
        ntsm_file="~{ntsm_output}"
        output_file="~{file_name}.summary.tsv"
        temp_file="~{file_name}.temp.tsv"

        # Extract values from the TSV file using awk
        total_reads=$(awk -F'\t' 'NR==2 {print $3}' ${readstat_file})
        total_bp=$(awk -F'\t' 'NR==3 {print $3}' ${readstat_file})
        total_Gbp=$(awk -F'\t' 'NR==4 {print $3}' ${readstat_file})
        total_min=$(awk -F'\t' 'NR==5 {print $3}' ${readstat_file})
        total_max=$(awk -F'\t' 'NR==6 {print $3}' ${readstat_file})
        total_mean=$(awk -F'\t' 'NR==7 {print $3}' ${readstat_file})
        quartile_25=$(awk -F'\t' 'NR==8 {print $3}' ${readstat_file})
        quartile_50=$(awk -F'\t' 'NR==9 {print $3}' ${readstat_file})
        quartile_75=$(awk -F'\t' 'NR==10 {print $3}' ${readstat_file})
        N25=$(awk -F'\t' 'NR==11 {print $3}' ${readstat_file})
        N50=$(awk -F'\t' 'NR==12 {print $3}' ${readstat_file})
        N75=$(awk -F'\t' 'NR==13 {print $3}' ${readstat_file})

        ntsm_values=$(awk -F'\t' '{print $(NF-1)"\t"$NF}' ${ntsm_file})
        ntsm_score=$(echo "$ntsm_values" | awk '{print $1}')
        ntsm_result=$(echo "$ntsm_values" | awk '{print $2}')
        
        # Create output file with a header
        echo -e "filename\ttotal_reads\ttotal_bp\ttotal_Gbp\tmin\tmax\tmean\tquartile_25\tquartile_50\tquartile_75\tN25\tN50\tN75\tntsm_score\tntsm_result" > ${output_file}
        echo -e "~{file_name}\t$total_reads\t$total_bp\t$total_Gbp\t$total_min\t$total_max\t$total_mean\t$quartile_25\t$quartile_50\t$quartile_75\t$N25\t$N50\t$N75\t$ntsm_score\t$ntsm_result" >> ${output_file}

        # Append methylation report if provided
        if [ "~{methylation_report}" != "/dev/null" ]; then
            paste ${output_file} ~{methylation_report} > ${temp_file}
            mv ${temp_file} ${output_file}
        fi
        
        # Log the generated summary file
        echo "Generated summary file: ${output_file}"
    >>>

    output {
        File summary_file = "~{file_name}.summary.tsv"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        docker: "ubuntu:latest"
        preemptible: preempts
    }
}
