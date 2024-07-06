version 1.0

workflow methylationWorkflow {
    input {
        Array[File] bam_files
        String sample_id
        Int addldisk = 20  # Default value set to 20 GB
    }

    Int read_size = ceil(size(bam_files, "GB"))
    Int final_disk_size = 3 * read_size + addldisk

    call runMethylationCheck {
        input:
            bam_files = bam_files,
            sample_id = sample_id,
            final_disk_size = final_disk_size,
            cpu = 4,
            memory_gb = 16
    }

    output {
        Array[File] methylation_outputs = runMethylationCheck.methylation_outputs
    }
}

task runMethylationCheck {
    input {
        Array[File] bam_files
        String sample_id
        Int final_disk_size
        Int cpu = 4
        Int memory_gb = 16
    }

    command <<<
        #!/bin/bash

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

          # Check for tags in the first 100 alignments
          ALIGNMENTS=$(samtools view "$BAM_FILE" 2>/dev/null | head -n 100)

          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}MM:Z:C\{0,10\}' > /dev/null; then
            MM="True"
          fi
          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}ML:B:C\{0,10\}' > /dev/null; then
            ML="True"
          fi
          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}fi:B:C\{0,10\}' > /dev/null; then
            FI="True"
          fi
          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}ri:B:C\{0,10\}' > /dev/null; then
            RI="True"
          fi
          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}fp:B:C\{0,10\}' > /dev/null; then
            FP="True"
          fi
          if echo "$ALIGNMENTS" | grep -o -m 1 '.\{0,10\}rp:B:C\{0,10\}' > /dev/null; then
            RP="True"
          fi

          # Output summary file
          SUMMARY_FILE="${FILE_NAME}_methylation_summary.txt"
          echo -e "$FILE_NAME\t$MM\t$ML\t$ALL_KINETICS_RESULT\t$KEEP_KINETICS_RESULT\t$HIFI_KINETICS_RESULT\t$PP_PRIMROSE\t$FI\t$RI\t$FP\t$RP" > "$SUMMARY_FILE"
        }

        # Process each BAM file
        for BAM_FILE in ~{sep=' ' bam_files}; do
          process_bam_file "$BAM_FILE"
        done

        echo "Processing completed"
    >>>

    runtime {
        docker: "quay.io/biocontainers/samtools:1.18--h50ea8bc_0"
        disks: "~{final_disk_size} GB"
        cpu: cpu
        memory: memory_gb + " GB"
    }

    output {
        Array[File] methylation_outputs = glob("*_methylation_summary.txt")
    }
}
