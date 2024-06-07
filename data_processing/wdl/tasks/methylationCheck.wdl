version 1.0

task runMethylationCheck {
    input {
        Array[File] bam_files
        String sample_id
    }


     command <<<
        set -e
        echo "Starting the task" || exit 1
        for bam_file in "~{sep=' ' bam_files}"; do
            echo "Processing $bam_file" || exit 1
            bam_basename=$(basename $bam_file .bam)
            methylation_output_file="${bam_basename}_methylation_summary.txt"
            if samtools view $bam_file 2>/dev/null | grep -m 1 "MM:Z:C+m" > /dev/null; then
                echo "$bam_file: yes" > $methylation_output_file
            else
                echo "$bam_file: no" > $methylation_output_file
            fi || exit 1
        done
    >>>

    runtime {
        docker: "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
    }

    output {
        Array[File] methylation_outputs = glob("*_methylation_summary.txt")
    }
}

workflow methylationWorkflow {
    input {
        Array[File] bam_files
        String sample_id
    }

    call runMethylationCheck {
        input:
            bam_files = bam_files,
            sample_id = sample_id
    }

    output {
        Array[File] methylation_outputs = runMethylationCheck.methylation_outputs
    }
}