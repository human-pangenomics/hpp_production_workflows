version 1.0

workflow runNucFreq {

    input {
        File input_bam
        File input_bam_bai
        File regions_bed
        String assembly_id
    }

    meta {
        description: "Calls NucFreq to identify regions in assembly with unexpected heterozygosity. See [NucFreq's documentation](https://github.com/mrvollger/NucFreq)"
    }
    
    parameter_meta {
        input_bam: "HiFi or ONT reads for a sample aligned to the sample's diploid assembly with Winnowmap or minimap2"
        input_bam_bai: "Index file for the input_bam"
        regions_bed: "Bed file of regions to investigate with NucFreq"
        assembly_id: "Name of assembly (for output file naming)"
    }


    call filter_bam {
        input:
            input_bam      = input_bam,
            input_bam_bai  = input_bam_bai,
            regions_bed    = regions_bed
    }

    call nucfreq {
        input:
            input_bam     = filter_bam.nucfreq_filt_bam,
            input_bam_bai = filter_bam.nucfreq_filt_bam_bai,
            regions_bed   = regions_bed,
            assembly_id   = assembly_id
    }

    call filter_nucfreq {
        input:
            variant_loci_bed = nucfreq.variant_loci_bed,
            assembly_id   = assembly_id
    }

  
    output {
        File nucplot_image    = nucfreq.nucplot_png
        File variant_loci_bed = nucfreq.variant_loci_bed
        File error_loci_bed   = filter_nucfreq.variant_clusters_bed
    }
}

task filter_bam {
    input{
        File input_bam
        File input_bam_bai
        File regions_bed 
        String sam_omit_flag = "2308"

        Int threadCount    = 8    
        Int memSizeGB      = 48
        Int addldisk       = 64    
        String dockerImage = "quay.io/biocontainers/samtools@sha256:9cd15e719101ae8808e4c3f152cca2bf06f9e1ad8551ed43c1e626cb6afdaa02" # 1.19.2--h50ea8bc_1
    }
    
    String file_prefix = basename(input_bam, ".bam")

    Int bam_size = ceil(size(input_bam, "GB"))
    Int final_disk_dize = 2*bam_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        samtools view \
            -F ~{sam_omit_flag}\
            --bam \
            --with-header \
            --regions-file ~{regions_bed} \
            --threads ~{threadCount} \
            -X ~{input_bam} ~{input_bam_bai} \
            -o ~{file_prefix}_nucfreq_filtered.bam

        samtools index \
            --threads ~{threadCount} \
            ~{file_prefix}_nucfreq_filtered.bam
  >>>  

  output {
    File nucfreq_filt_bam     = "~{file_prefix}_nucfreq_filtered.bam"
    File nucfreq_filt_bam_bai = "~{file_prefix}_nucfreq_filtered.bam.bai"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}


task nucfreq {
    input{
        File input_bam
        File input_bam_bai
        File regions_bed 
        String assembly_id 

        String sam_omit_flag = "2308"
        String otherArgs   = ""

        Int threadCount    = 20   
        Int memSizeGB      = 80
        Int addldisk       = 64    
        String dockerImage = "humanpangenomics/nucfreq@sha256:6f2f981892567f2a8ba52ba20e87f98e6ca770ea3f4d5430bf67a26673c8f176" 
    }

    Int bam_size = ceil(size(input_bam, "GB"))
    Int final_disk_dize = bam_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        ## soft link bam and bai to cwd so they are in the same directory
        ln -s ~{input_bam} input.bam
        ln -s ~{input_bam_bai} input.bam.bai

        ## run nucfreq: find loci with heterozygosity; create plot of regions
        python /opt/nucfreq/NucPlot.py \
            -t ~{threadCount} \
            --bed ~{regions_bed} \
            --obed ~{assembly_id}_variant_loci.bed \
            input.bam \
            ~{assembly_id}_nucplot \
            ~{otherArgs}
  >>>  

  output {
    File variant_loci_bed = "~{assembly_id}_variant_loci.bed"
    File nucplot_png      = "~{assembly_id}_nucplot.png"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}


task filter_nucfreq {
    input{
        File variant_loci_bed
        String assembly_id 

        String otherArgs   = ""

        Int threadCount    = 4   
        Int memSizeGB      = 16
        Int addldisk       = 32    
        String dockerImage = "rocker/verse@sha256:56e60da5b006e1406967e58ad501daaba567d6836029aee94ed16ba1965554f0" # 4.3.1
    }

    Int bed_size = ceil(size(variant_loci_bed, "GB"))
    Int final_disk_dize = bed_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        wget https://raw.githubusercontent.com/emics57/nucfreqPipeline/89b1b9a980f78b63e2a79ca2e122269bf284df41/nucfreq_filtering_migalab.R

        Rscript nucfreq_filtering_migalab.R \
            ~{variant_loci_bed} \
            ~{assembly_id}_nucfreq_errors.bed \
            ~{otherArgs}

  >>>  

  output {
    File variant_clusters_bed = "~{assembly_id}_nucfreq_errors.bed"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
