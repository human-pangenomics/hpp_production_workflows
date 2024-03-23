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

    call nucfreq as nucfreq_vanilla {
        input:
            input_bam     = filter_bam.nucfreq_filt_bam,
            input_bam_bai = filter_bam.nucfreq_filt_bam_bai,
            regions_bed   = regions_bed,
            assembly_id   = assembly_id,
            tag           = "discordant"
    }

    call filter_nucfreq {
        input:
            nucfreq_loci_bed = nucfreq_vanilla.nucfreq_loci_bed,
            assembly_id   = assembly_id
    }

    call nucfreq as nucfreq_all_loci_in_regions {
        input:
            input_bam     = filter_bam.nucfreq_filt_bam,
            input_bam_bai = filter_bam.nucfreq_filt_bam_bai,
            regions_bed   = regions_bed,
            assembly_id   = assembly_id,
            otherArgs     = "--minobed 0",
            tag           = "all"
    }


    output {
        File nucplot_image_tar  = nucfreq_vanilla.nucplot_images
        File variant_loci_bed   = nucfreq_vanilla.nucfreq_loci_bed
        File all_loci_bed       = nucfreq_all_loci_in_regions.nucfreq_loci_bed
        File error_clusters_bed = filter_nucfreq.variant_clusters_bed        
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

        String tag = ""
        String otherArgs   = ""

        Int threadCount    = 4   
        Int memSizeGB      = 32
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

        # Split bed to one file per row. Run on one row at a time
        # Create a directory to store split BED files
        mkdir -p split_beds
        mkdir -p split_beds_out
        mkdir -p output_plots


        ## run nucfreq: find loci with heterozygosity; create plots for each region
        while IFS=$'\t' read -r chrom start end rest; do
            
            FILE_NAME="${chrom}_${start}_${end}.bed"
            echo -e "$chrom\t$start\t$end\t$rest" > "split_beds/$FILE_NAME"

            python /opt/nucfreq/NucPlot.py \
                -t ~{threadCount} \
                --bed "split_beds/$FILE_NAME" \
                --obed "split_beds_out/$FILE_NAME" \
                input.bam \
                "output_plots/~{assembly_id}_${chrom}_${start}_${end}.png" \
                ~{otherArgs}

        done < ~{regions_bed}

        
        # Process the first file fully, including the header
        head -n 1 split_beds_out/$(ls split_beds_out | head -n 1) > ~{assembly_id}_nucfreq_loci_bed

        # Concatenate the rest of the files without the header and then sort
        for file in split_beds_out/*.bed; do
            tail -n +2 "$file"
        done | sort -k1,1 -k2,2n >> "~{assembly_id}_nucfreq_loci_~{tag}.bed"

        ## tar.gz individual plots 
        tar -czvf "~{assembly_id}_nucfreq_plots.tar.gz" output_plots

  >>>  

  output {
    File nucfreq_loci_bed = "~{assembly_id}_nucfreq_loci_~{tag}.bed"
    File nucplot_images   = "~{assembly_id}_nucfreq_plots.tar.gz"
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
        File nucfreq_loci_bed
        String assembly_id 

        String otherArgs   = ""

        Int threadCount    = 4   
        Int memSizeGB      = 16
        Int addldisk       = 32    
        String dockerImage = "rocker/verse@sha256:56e60da5b006e1406967e58ad501daaba567d6836029aee94ed16ba1965554f0" # 4.3.1
    }

    Int bed_size = ceil(size(nucfreq_loci_bed, "GB"))
    Int final_disk_dize = bed_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        wget https://raw.githubusercontent.com/emics57/nucfreqPipeline/21b3395a7f285962aae9e881db2514e03601c5db/nucfreq_filtering_migalab.R

        Rscript nucfreq_filtering_migalab.R \
            ~{nucfreq_loci_bed} \
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
