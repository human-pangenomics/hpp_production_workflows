version 1.0

workflow rm2bed_workflow {

    call rm2bed 

    output {
        File rm_bed = rm2bed.rm_bed
    }
}



task rm2bed {

    input {
        String sample_name
        String output_file_tag
        File rm_out_file

        Int memSizeGB = 4
        Int diskSizeGB = 32
    }

    String output_bed_fn = "${sample_name}.${output_file_tag}_rm.bed"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Call RM2Bed
        RM2Bed.py ~{rm_out_file}

        ## sort output bed file
        sort -V -k1,1 -k2,2 *_rm.bed ~{output_bed_fn}

    >>>

    output {

        File rm_bed  = output_bed_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/rm2bed@sha256:ce87d2ff57f4a1755b79a0c1aee3fd0be47a18918f04c257d509305db6ab82e8"
        preemptible: 1
    }
}