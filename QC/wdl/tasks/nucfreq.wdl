version 1.0

workflow runNucFreq {

    input {
        File input_bam
        File input_bam_bai
        File regions_bed
        File assembly_fasta
    }

    meta {
        description: "Calls NucFreq to identify regions in assembly with unexpected heterozygosity. See [NucFreq's documentation](https://github.com/mrvollger/NucFreq)"
    }
    
    parameter_meta {
        input_bam: "HiFi or ONT reads for a sample aligned to the sample's diploid assembly with Winnowmap or minimap2"
        input_bam_bai: "Index file for the input_bam"
        assembly_fasta: "Assembly that reads were aligned against. Can be fasta, fasta.gz, fa, or ga.gz. Used for creating genome regions to split."
        regions_bed: "Bed file of regions in which to output NucFreq plots"
    }

    String output_prefix = basename(input_bam, ".bam")

    call filter_bam {
        input:
            input_bam      = input_bam,
            input_bam_bai  = input_bam_bai
    }

    ## Call nucreq in just regions_bed to get nucfreq plots
    call nucfreq {
        input:
            input_bam     = filter_bam.nucfreq_filt_bam,
            input_bam_bai = filter_bam.nucfreq_filt_bam_bai,
            regions_bed   = regions_bed
    } 


    call create_genome_beds {
        input:
            inputFasta = assembly_fasta
    }

    scatter (genome_bed in create_genome_beds.genome_beds) {
        call nucfreq_counts {
            input:
                input_bam     = filter_bam.nucfreq_filt_bam,
                input_bam_bai = filter_bam.nucfreq_filt_bam_bai,
                regions_bed   = genome_bed
        }

        call create_nucfreq_output {
            input:
                nucfreq_counts_txt = nucfreq_counts.nucfreq_counts_txt
        }

    }
    
    
    call combine_beds as combine_first_bedgraph {
        input:
            beds         = create_nucfreq_output.nucfreq_first_bedgraph,
            prefix       = output_prefix,
            tag          = "nucfreq_first"
    }

    call combine_beds as combine_second_bedgraph {
        input:
            beds         = create_nucfreq_output.nucfreq_second_bedgraph,
            prefix       = output_prefix,
            tag          = "nucfreq_second"            
    }

    call combine_beds as combine_nucfreq_bed {
        input:
            beds         = create_nucfreq_output.nucfreq_bed,
            prefix       = output_prefix,
            tag          = "nucfreq"            
    }    


    call filter_nucfreq {
        input:
            nucfreq_loci_bed = combine_nucfreq_bed.output_bed
    }

    call bedgraph_to_bigwig as first_allele_bw {
        input:
            bedgraph    = combine_first_bedgraph.output_bed,
            chrom_sizes = create_genome_beds.chrom_sizes
    }
     
    call bedgraph_to_bigwig as second_allele_bw {
        input:
            bedgraph    = combine_second_bedgraph.output_bed,
            chrom_sizes = create_genome_beds.chrom_sizes
    }       

    output {
        ## If regions were passed as inputs
        File nucplot_image_tar        = nucfreq.nucplot_images
        File nucfreq_all_bed          = combine_nucfreq_bed.output_bed       
        File error_clusters_bed       = filter_nucfreq.variant_clusters_bed        
        File first_allele_bigwig      = first_allele_bw.bigwig
        File second_allele_bigwig     = second_allele_bw.bigwig
    }
}

task filter_bam {
    input{
        File input_bam
        File input_bam_bai
        File? regions_bed 
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

        if [ ! -f "~{regions_bed}" ]
        then
            REGIONS_ARG=""
        else
            REGIONS_ARG="--regions-file ~{regions_bed}"
        fi

        samtools view \
            -F ~{sam_omit_flag}\
            --bam \
            --with-header \
            $REGIONS_ARG \
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

        String tag = ""
        String otherArgs   = ""

        Int threadCount    = 4   
        Int memSizeGB      = 32
        Int addldisk       = 64    
        String dockerImage = "humanpangenomics/nucfreq@sha256:6f2f981892567f2a8ba52ba20e87f98e6ca770ea3f4d5430bf67a26673c8f176" 
    }

    String file_prefix = basename(input_bam, ".bam")

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
                "output_plots/~{file_prefix}_${chrom}_${start}_${end}.png" \
                ~{otherArgs}

        done < ~{regions_bed}

        
        # Process the first file fully, including the header
        head -n 1 split_beds_out/$(ls split_beds_out | head -n 1) > "~{file_prefix}_nucfreq_loci.bed"

        # Concatenate the rest of the files without the header and then sort
        for file in split_beds_out/*.bed; do
            tail -n +2 "$file"
        done | sort -k1,1 -k2,2n >> "~{file_prefix}_nucfreq_regions_loci.bed"

        ## tar.gz individual plots 
        tar -czvf "~{file_prefix}_nucfreq_plots.tar.gz" output_plots

  >>>  

  output {
    File nucfreq_loci_bed = "~{file_prefix}_nucfreq_regions_loci.bed"
    File nucplot_images   = "~{file_prefix}_nucfreq_plots.tar.gz"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}

task define_genome_split_beds {
    input{
        File inputFasta
        
        Int region_size    = 50000000

        Int threadCount    = 2   
        Int memSizeGB      = 16
        Int diskSize       = 32
        String dockerImage = "quay.io/biocontainers/samtools@sha256:9cd15e719101ae8808e4c3f152cca2bf06f9e1ad8551ed43c1e626cb6afdaa02" # 1.19.2--h50ea8bc_1
    }

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        inputFastaFN=$(basename -- "~{inputFasta}")

        ## first check if inputFasta needs to be unzipped
        if [[ $inputFastaFN =~ \.gz$ ]]; then
            cp ~{inputFasta} .
            gunzip -f $inputFastaFN
            inputFastaFN="${inputFastaFN%.gz}"
        else
            ln -s ~{inputFasta}
        fi 

        ## get contig/scaffold sizes from genome assembly
        samtools faidx "$inputFastaFN" 
        cut -f1,2 "${inputFastaFN}.fai" > sizes.genome

        mkdir genomic_windows

        ## split genome into regions; output bed file for each region (to scatter on)
        while IFS=$'\t' read -r CHROM SIZE; do
        
            # find number of windows. Add 1 if there is a remainder.
            NUMWINS=$((SIZE / ~{region_size} + (SIZE % ~{region_size} > 0 ? 1 : 0)))

            for ((i=1; i<=NUMWINS; i++)); do
                START=$(( (i - 1) * ~{region_size} + 1 ))
                END=$(( i * ~{region_size} ))
                END=$(( END < SIZE ? END : SIZE ))

                FILE_NAME="${CHROM}_${START}_${END}.bed"
                echo -e "$CHROM\t$START\t$END\tchrom_split" > "genomic_windows/$FILE_NAME"
            done
        done < sizes.genome


  >>>  

  output {
    Array[File] genomic_windows_beds = glob("genomic_windows/*.bed")
    File chrom_sizes = "sizes.genome"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}

task create_genome_beds {
    input{
        File inputFasta

        Int threadCount    = 2   
        Int memSizeGB      = 16
        Int diskSize       = 32
        String dockerImage = "quay.io/biocontainers/samtools@sha256:9cd15e719101ae8808e4c3f152cca2bf06f9e1ad8551ed43c1e626cb6afdaa02" # 1.19.2--h50ea8bc_1
    }

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        inputFastaFN=$(basename -- "~{inputFasta}")

        ## first check if inputFasta needs to be unzipped
        if [[ $inputFastaFN =~ \.gz$ ]]; then
            cp ~{inputFasta} .
            gunzip -f $inputFastaFN
            inputFastaFN="${inputFastaFN%.gz}"
        else
            ln -s ~{inputFasta}
        fi 

        ## get contig/scaffold sizes from genome assembly
        samtools faidx "$inputFastaFN" 
        cut -f1,2 "${inputFastaFN}.fai" > unsorted_sizes.genome

        sort -k1,1 unsorted_sizes.genome > sizes.genome

        ## create bed files with one chromosome per file (for scattering on)
        mkdir -p beds

        ## name with original line number to make aggregating easier
        line_number=1

        # Read through each line of the file
        while IFS= read -r line; do

            chrom=$(echo "$line" | awk '{print $1}')
            size=$(echo "$line" | awk '{print $2}')
        
            echo -e "${chrom}\t0\t${size}" > "beds/genome_region_${line_number}.bed"

            ((line_number++))

        done < sizes.genome          
  >>>  

  output {
    Array[File] genome_beds  = glob("beds/*.bed")
    File chrom_sizes         = "sizes.genome"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}

task nucfreq_counts {
    input{
        File input_bam
        File input_bam_bai
        File regions_bed

        Int threadCount    = 8   
        Int memSizeGB      = 32
        Int addldisk       = 128   
        String dockerImage = "quay.io/biocontainers/rustybam@sha256:0c31acc94fe676fd7d853da74660187d9a146acbacb4266abd2ec559fd5641a3" # 0.1.33--h756b843_0
    }
    String file_prefix = basename(regions_bed, ".bed")

    Int bam_size = ceil(size(input_bam, "GB"))
    Int final_disk_dize = bam_size * 2 + addldisk

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        
        ## soft link bam and bai to cwd so they are in the same directory
        ln -s ~{input_bam} input.bam
        cp ~{input_bam_bai} input.bam.bai

        rustybam nucfreq \
            --bed ~{regions_bed} \
            input.bam \
            > "~{file_prefix}_nucfreq_counts.txt" 2> stderr.txt 

        ## remove header line        
        sed -i '1d' "~{file_prefix}_nucfreq_counts.txt"

  >>>

  output {
    File nucfreq_counts_txt  = "~{file_prefix}_nucfreq_counts.txt"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}


task create_nucfreq_output {
    input{
        File nucfreq_counts_txt

        Int threadCount    = 4   
        Int memSizeGB      = 32
        Int addldisk       = 64    
        String dockerImage = "quay.io/biocontainers/mawk@sha256:793a6f7e2f641288585b0714d097c5c6b146d5d54c72eb000de2dc7bef24a784" # 1.3.4--h031d066_7
    }

    String file_prefix = basename(nucfreq_counts_txt, ".txt")

    Int file_size = ceil(size(nucfreq_counts_txt, "GB"))
    Int final_disk_dize = file_size * 3 + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        out_name="~{file_prefix}"

        mawk '{
            max = $4; second_max = 0;

            for (i = 5; i <= 7; i++) {
                if ($i > max) { second_max = max; max = $i; }
                else if ($i > second_max) { second_max = $i; }
            }

            # Construct file names
            first_file = "'"${out_name}"'_nucfreq_first.bedGraph";
            second_file = "'"${out_name}"'_nucfreq_second.bedGraph";
            combined_file = "'"${out_name}"'_nucfreq.bed";

            print $1, $2, $3, max >> first_file;
            print $1, $2, $3, second_max >> second_file;
            print $1, $2, $3, max, second_max >> combined_file;
        }

        END {
            close(first_file);
            close(second_file);
            close(combined_file);

        }' < ~{nucfreq_counts_txt}

        
  >>>

  output {
    File nucfreq_first_bedgraph  = "~{file_prefix}_nucfreq_first.bedGraph"
    File nucfreq_second_bedgraph = "~{file_prefix}_nucfreq_second.bedGraph"
    File nucfreq_bed             = "~{file_prefix}_nucfreq.bed"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}

task combine_beds {
    input {
        Array[File] beds        
        String prefix = "sample"
        String tag    = "combined"

        Int memSizeGB   = 64
        Int threadCount = 4
        Int addldisk = 128
    }

    Int file_size = ceil(size(beds, "GB"))
    Int final_disk_dize = file_size * 2 + addldisk

    String out_bed_fn  = "~{prefix}_~{tag}.bed"

    command <<<
        set -eux -o pipefail

        ## Combine scattered results into one file
        LC_ALL=C

        cat ~{sep=" " beds} | sort -k1,1 -k2,2n > ~{out_bed_fn}

    >>>

    output {
        File output_bed  = out_bed_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "juklucas/alphasat_summarize@sha256:872277a3a780c676d32f61045043b638e07889e63faa0de054e5f3d7d362f7b4"
        preemptible: 1
    }

}

task filter_nucfreq {
    input{
        File nucfreq_loci_bed

        String otherArgs   = ""

        Int threadCount    = 4   
        Int memSizeGB      = 16
        Int addldisk       = 32    
        String dockerImage = "rocker/verse@sha256:56e60da5b006e1406967e58ad501daaba567d6836029aee94ed16ba1965554f0" # 4.3.1
    }
    String file_prefix = basename(nucfreq_loci_bed, ".bed")

    Int bed_size = ceil(size(nucfreq_loci_bed, "GB"))
    Int final_disk_dize = bed_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        wget https://raw.githubusercontent.com/emics57/nucfreqPipeline/21b3395a7f285962aae9e881db2514e03601c5db/nucfreq_filtering_migalab.R

        Rscript nucfreq_filtering_migalab.R \
            ~{nucfreq_loci_bed} \
            ~{file_prefix}_errors.bed \
            ~{otherArgs}
  >>>  

  output {
    File variant_clusters_bed = "~{file_prefix}_errors.bed"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}

task bedgraph_to_bigwig {
    input{
        File bedgraph
        File chrom_sizes

        Int threadCount    = 4    
        Int memSizeGB      = 12
        Int addldisk       = 64    
        String dockerImage = "quay.io/biocontainers/ucsc-bedgraphtobigwig@sha256:9a5a150acf6af3910d939396e928dc3d9468d974624eef7fc74ab6e450c12466" # 455--h2a80c09_1
    }
    
    String file_prefix = basename(bedgraph, ".bed")

    Int bed_size = ceil(size(bedgraph, "GB"))
    Int final_disk_dize = bed_size + addldisk

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        bedGraphToBigWig ~{bedgraph} ~{chrom_sizes} ~{file_prefix}.bw
  >>>  

  output {
    File bigwig     = "~{file_prefix}.bw"
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + final_disk_dize + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
