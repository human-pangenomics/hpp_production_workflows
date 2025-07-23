version 1.0

workflow run_primrose {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Outputs hifi bam file with 5-Methylcytosine (5mC) predictions of each CpG in PacBio HiFi reads with [primrose](https://github.com/PacificBiosciences/primrose)"
    }

    input {
        Array[File] input_bams
    }

    scatter (input_bam in input_bams) {
        call primrose {
            input:
                input_bam = input_bam
        }
    }

    output {
        Array[File] output_5mc_bam = primrose.output_bam
    }  
}

task primrose {
    input {
        File input_bam
        String? output_name
        
        Int memSizeGB   = 64
        Int threadCount = 8
        Int addlDisk    = 50        
        Int preempts    = 3
        String extraArgs = "--keep-kinetics"
    }

    parameter_meta {
        input_bam: "Input for primrose are PacBio HiFi reads with kinetics. File should be named *.hifi_reads.bam"
        output_name: "Optional custom output filename. If not provided, will use input basename with .5mc.hifi_reads.bam suffix"
    }

    # Generate default output name if not provided
    String input_basename = basename(input_bam, ".bam")
    String default_output_name = if (sub(input_basename, ".*\\.hifi_reads$", "MATCH") == "MATCH")
        then sub(input_basename, "\\.hifi_reads$", ".5mc.hifi_reads.bam")
        else input_basename + ".5mc.hifi_reads.bam"
    
    String final_output_name = select_first([output_name, default_output_name])

    # Estimate disk size required
    Int input_bam_size   = ceil(size(input_bam, "GB"))    
    Int final_disk_size  = input_bam_size * 2 + addlDisk

    command <<<
        set -eux -o pipefail

        ## Pass extra arguments if extraArgs is set, if not just pass empty string
        if [ -z "~{extraArgs}" ]
        then
            EXTRA_ARGS=""
        else
            EXTRA_ARGS="~{extraArgs}"
        fi
            
        primrose ${EXTRA_ARGS} ~{input_bam} ~{final_output_name}
    >>>

    output {
        File output_bam = final_output_name
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "humanpangenomics/primrose@sha256:93ed60851f6a43e46e20b4075b0fea146e8ecc0b86d00b6800b987f0253fcd2e"
        preemptible: preempts
    }
}