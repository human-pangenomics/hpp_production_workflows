version 1.0

import "primrose.wdl" as primrose_wdl

workflow run_primrose_with_kinetics_fix {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        maintainer: "Andrew Blair"
        description: "Fixes kinetics arrays in HiFi BAM files (for files demultiplexed with lima <2.5) then outputs hifi bam file with 5-Methylcytosine (5mC) predictions of each CpG in PacBio HiFi reads with [primrose](https://github.com/PacificBiosciences/primrose). Only run if input bam has kinetics tags which are incorrectly trimmed by lima."
    }

    input {
        File input_bam
        Boolean drop_mm = true
        String? output_name
        String sample_id  # Only used to help name outputs if output_name is not given
    }

    # Use sample_id only if output_name wasn't explicitly provided
    String default_output_name = sample_id + ".fixed_kinetics.bam"
    String final_output_name = select_first([output_name, default_output_name])

    call trim_hifi_barcode_kinetics {
        input:
            input_bam = input_bam,
            drop_mm = drop_mm
    }

    call primrose_wdl.primrose as primrose_task {
        input:
            input_bam = trim_hifi_barcode_kinetics.output_bam,
            output_name = final_output_name
    }

    output {
        File output_5mc_bam = primrose_task.output_bam
    }

    parameter_meta {
        input_bam: "Input PacBio HiFi reads BAM file with barcode and kinetics tags"
        drop_mm: "Whether to drop MM/ML methylation tags (default: true)"
        output_name: "Optional custom name for the final primrose BAM output"
        sample_id: "Sample identifier used for default output BAM naming"
    }
}

task trim_hifi_barcode_kinetics {
    input {
        File input_bam
        Boolean drop_mm = true

        Int memSizeGB   = 32
        Int threadCount = 4
        Int addlDisk    = 100        
        Int preempts    = 1
    }

    parameter_meta {
        input_bam: "Input PacBio HiFi reads with kinetics that were demultiplexed with lima <2.5. Should have bx tags and kinetics arrays that need trimming."
        drop_mm: "Remove existing MM/ML methylation tags"
    }

    Int input_bam_size   = ceil(size(input_bam, "GB"))    
    Int final_disk_size  = input_bam_size * 2 + addlDisk

    command <<<
        set -eux -o pipefail

        # Extract base filename without .bam extension
        FILENAME=$(basename -- "~{input_bam}" .bam)
        
        # Set up drop_mm flag
        if [ "~{drop_mm}" == "true" ]
        then
            DROP_MM_FLAG="--drop_mm"
        else
            DROP_MM_FLAG=""
        fi

        # Run the kinetics trimming script
        python /opt/trim_hifi_barcode_kinetics.py \
            ${DROP_MM_FLAG} \
            --threads ~{threadCount} \
            ~{input_bam} \
            ${FILENAME}.kinetics_trimmed.bam
    >>>

    output {
        File output_bam = glob("*.kinetics_trimmed.bam")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "humanpangenomics/trim_hifi_barcode_kinetics:1.0.0--8f4990344577ac221ceb0d546f61fc305d7e785b"
        preemptible: preempts
    }
}
