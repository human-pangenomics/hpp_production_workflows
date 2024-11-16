version 1.0

workflow bed_to_bigbed_wf {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Convert bed to bigbed format for browser viewing."
    }

    parameter_meta {
        input_bed: "Input bed to convert to bigbed format."
        auto_sql: "autoSql (*.as) file which specifies bed column definitions. See https://genomewiki.ucsc.edu/index.php/AutoSql"
        assembly_fa: "Assembly fasta (can be gzipped) to extract chromosome sizes from."
        type_str: "bedN[+[P]] where N is number of fields and P is the number of bedPlus fields. For example bed6 or bed6+3"
        tab_delim: "If bed file has whitespace (not tab) delimeters, enter false"
    }

    input {
        File input_bed
        File auto_sql
        File assembly_fa
        String type_str   = "bed6"
        Boolean tab_delim = true 
    }

    call create_chrom_sizes {
        input:
            assembly_fa = assembly_fa
    }

    call bed_to_bigbed {
        input:
            input_bed   = input_bed,
            auto_sql    = auto_sql,
            type_str    = type_str,
            tab_delim   = tab_delim,
            chrom_sizes = create_chrom_sizes.chrom_sizes
    }

    output {
        File output_bigbed = bed_to_bigbed.output_bb
    }
}

task create_chrom_sizes {
    input{
        File assembly_fa

        Int thread_count    = 2   
        Int mem_size_gb     = 16
        Int disk_size       = 32
        String docker_image = "quay.io/biocontainers/samtools@sha256:9cd15e719101ae8808e4c3f152cca2bf06f9e1ad8551ed43c1e626cb6afdaa02" # 1.19.2--h50ea8bc_1
        Int preemptible     = 1
    }

    String output_prefix = sub(basename(assembly_fa), "\.(fa|fasta)(\.gz)?$", "")

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        assembly_fn=$(basename -- "~{assembly_fa}")

        ## first check if assembly_fa needs to be unzipped
        if [[ $assembly_fn =~ \.gz$ ]]; then
            cp ~{assembly_fa} .
            gunzip -f $assembly_fn
            assembly_fn="${assembly_fn%.gz}"
        else
            ln -s ~{assembly_fa}
        fi 

        ## get contig/scaffold sizes from genome assembly
        samtools faidx "$assembly_fn" 
        cut -f1,2 "${assembly_fn}.fai" > "~{output_prefix}.chrom_sizes"
    >>>  

    output {
        File chrom_sizes = "~{output_prefix}.chrom_sizes"
    }

    runtime {
        cpu: thread_count
        memory: mem_size_gb + " GB"    
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}

task bed_to_bigbed {
    input {
        File input_bed
        File auto_sql
        Boolean tab_delim
        String type_str
        File chrom_sizes
        
        Int thread_count    = 2   
        Int mem_size_gb     = 16
        Int disk_size       = 32
        String docker_image = "quay.io/biocontainers/ucsc-bedtobigbed@sha256:be939493b8ed5f504a5a2f1ac1bcbcf6c30a6d7378b482f65bba3be596035e9e"
        Int preemptible     = 1
    }
    
    String delimeter_arg = if tab_delim then "-tab" else ""
    String output_prefix = basename(input_bed, ".bed")

    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        bedToBigBed \
            ~{delimeter_arg} \
            -type=~{type_str} \
            -as=~{auto_sql} \
            ~{input_bed} \
            ~{chrom_sizes} \
            ~{output_prefix}.bb
    >>>  

    output {
        File output_bb = "~{output_prefix}.bb"
    }

    runtime {
        cpu: thread_count
        memory: mem_size_gb + " GB"    
        disks: "local-disk " + disk_size + " SSD"
        docker: docker_image
        preemptible: preemptible
    }
}