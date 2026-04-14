version 1.0

import "../tasks/ntsm.wdl" as ntsm_tasks

workflow ntsm_eval_wf {
    input {
        Array[File] count_files
        String output_prefix
    }

    parameter_meta {
        count_files: "All ntsm count files to evaluate together (from ntsm_count task across all samples and read types)."
        output_prefix: "Prefix for the output eval TSV (e.g. batch name or project name)."
    }

    call ntsm_tasks.ntsm_eval {
        input:
            count_files   = count_files,
            output_prefix = output_prefix
    }

    output {
        File ntsm_eval_out = ntsm_eval.ntsm_eval_out
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Runs ntsmEval across all ntsm count files from a batch of QC runs to detect sample swaps."
    }
}
