version 1.0

import "../tasks/read_stats.wdl" as read_stats_wf
import "../tasks/ntsm.wdl" as ntsm_check

workflow hifi_qc_wf {
    input {
        Array[File] hifi_reads
        Array[File] other_reads
        String sample_id
    }

    parameter_meta {
        hifi_reads: "Reads to QC. Can be fastq, fastq.gz, or bam."
        other_reads: "Trusted reads to compare against to check for sample swaps. Can be fastq, fastq.gz, or bam."
        sample_id: "Id to use for output naming."
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
    }
    
    call read_stats_wf.runReadStats as read_stats {
        input: 
            reads          = hifi_reads,
            identifier     = sample_id
    }

    call ntsm_check.ntsm_workflow as ntsm_wf {
        input:
            input_reads_1  = hifi_reads,
            input_reads_2  = other_reads,
            sample_id      = sample_id,
            read_1_type    = "hifi",
            read_2_type    = "other"
    }


    output {
        ## readstats output
        File hifi_readstat_tar    = read_stats.ReadStatsTarball 
        File hifi_readstat_report = read_stats.ReadStatsReport

        ## ntsm output
        File hifi_ntsm_counts     = ntsm_wf.ntsm_count_1
        File ext_ntsm_counts      = ntsm_wf.ntsm_count_2
        File ntsm_eval            = ntsm_wf.ntsm_eval_out
    }
}

