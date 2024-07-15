version 1.0

workflow methylationWorkflow {
    input {
        Array[File] bam_files
        String sample_id
    }

    Int base_disk_size = 16  # Assuming a base disk size of 16 GB for temporary storage
    Int total_size = ceil(size(bam_files, "GB"))
    Int final_disk_size = total_size + base_disk_size

    call runMethylationCheck {
        input:
            bam_files = bam_files,
            sample_id = sample_id,
            memSizeGB = 4,
            threadCount = 1,
            disk_size = final_disk_size,
            preempts = 2
    }

    output {
        File methylation_summary = runMethylationCheck.methylation_summary
    }
}

task runMethylationCheck {
    input {
        Array[File] bam_files
        String sample_id
        Int memSizeGB
        Int threadCount
        Int disk_size
        Int preempts
    }

    command <<<
        #!/bin/bash

        set -eux

        # Define column headers
        headers="bam_file\tMM_tag\tML_tag\tall_kinetics_flag\tkeep_kinetics_flag\thifi_kinetics_tag\tPP_PRIMROSE\tfi_tag\tri_tag\tfp_tag\trp_tag"

        # Output summary file
        SUMMARY_FILE="~{sample_id}_methylation_summary.tsv"

        echo "Creating summary file: ${SUMMARY_FILE}"

        # Write headers to the summary file
        echo -e "$headers" > "$SUMMARY_FILE"
        echo "Headers written to summary file"

        # Function to process each BAM file
        process_bam_file() {
          local BAM_FILE=$1

          # Extract the last part of the path after the last '/'
          FILE_NAME=$(basename "$BAM_FILE")

          # Initialize tag variables
          MM="False"
          ML="False"
          ALL_KINETICS_RESULT="False"
          KEEP_KINETICS_RESULT="False"
          HIFI_KINETICS_RESULT="False"
          PP_PRIMROSE="False"
          FI="False"
          RI="False"
          FP="False"
          RP="False"

          # Check for "--all-kinetics" flag in the header
          if samtools view -H "$BAM_FILE" 2>/dev/null | grep -q -- "--all-kinetics"; then
            ALL_KINETICS_RESULT="True"
          fi

          # Check for "--keep-kinetics" flag in the header
          if samtools view -H "$BAM_FILE" 2>/dev/null | grep -q -- "--keep-kinetics"; then
            KEEP_KINETICS_RESULT="True"
          fi

          # Check for "--hifi-kinetics" flag in the header
          if samtools view -H "$BAM_FILE" 2>/dev/null | grep -q -- "--hifi-kinetics"; then
            HIFI_KINETICS_RESULT="True"
          fi

          # Check for "PP:primrose" in the header
          if samtools view -H "$BAM_FILE" 2>/dev/null | grep -q 'PP:primrose'; then
            PP_PRIMROSE="True"
          fi

          # Check for tags in the first 100 BAM_READS
          BAM_READS=$(samtools view "$BAM_FILE" 2>/dev/null | head -n 100)

          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}MM:Z:C\{0,10\}' > /dev/null; then
            MM="True"
          fi
          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}ML:B:C\{0,10\}' > /dev/null; then
            ML="True"
          fi
          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}fi:B:C\{0,10\}' > /dev/null; then
            FI="True"
          fi
          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}ri:B:C\{0,10\}' > /dev/null; then
            RI="True"
          fi
          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}fp:B:C\{0,10\}' > /dev/null; then
            FP="True"
          fi
          if echo "$BAM_READS" | grep -o -m 1 '.\{0,10\}rp:B:C\{0,10\}' > /dev/null; then
            RP="True"
          fi

          # Append results to the summary file
          echo -e "$FILE_NAME\t$MM\t$ML\t$ALL_KINETICS_RESULT\t$KEEP_KINETICS_RESULT\t$HIFI_KINETICS_RESULT\t$PP_PRIMROSE\t$FI\t$RI\t$FP\t$RP" >> "$SUMMARY_FILE"
          echo "Processed BAM file: $BAM_FILE"
        }

        # Process each BAM file
        for BAM_FILE in ~{sep=' ' bam_files}; do
          if [ -f "$BAM_FILE" ]; then
            process_bam_file "$BAM_FILE"
          else
            echo "Error: File $BAM_FILE does not exist or is not accessible." >&2
          fi
        done

        echo "Summary file created: $SUMMARY_FILE"
    >>>

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        docker: "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
        preemptible: preempts
    }

    output {
        File methylation_summary = "~{sample_id}_methylation_summary.tsv"
    }
}
