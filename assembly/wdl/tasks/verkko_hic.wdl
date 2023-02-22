version 1.0


workflow run_verkko_hic {

    call verkko_hic
    
    output{
        File hic_phased_run_folder = verkko_hic.hic_phased_run_folder
    }
}

task verkko_hic {
    
    input {
        File run_folder
        Array[File] input_hifi
        Array[File] input_nanopore
        Array[File] input_hic      
        String name = "assembly"
        String tag  = "verkko_hic"

        Int threadCount = 46
        Int memSizeGB   = 290
        Int diskSizeGB  = 2500
        Int preemptible = 0
    }

    command <<<

        set -eux -o pipefail

        ## localize nanopore reads to one directory
        ont_files=(~{sep=" " input_nanopore})

        mkdir -p inputs/ont 
        
        for ont_file in ${ont_files[@]};
        do
          cp $ont_file inputs/ont/
        done

        ## localize HiFi reads to one directory
        hifi_files=(~{sep=" " input_hifi})
        
        mkdir  -p inputs/hifi 

        for hifi_file in ${hifi_files[@]};
        do
            ## Have to copy so I can set timestamp
            cp $hifi_file inputs/hifi/
        done


        ## Set modtime to a few years ago so it doesn't trigger snakemake
        touch -a -m -t 202001011205.02 inputs/hifi/*
        touch -a -m -t 202001011205.02 inputs/ont/*


        ## localize HiC reads to one directory
        hic_files=(~{sep=" " input_hic})
        
        mkdir  -p inputs/hic 

        for hic_file in ${hic_files[@]};
        do
            ## Have to copy so I can set timestamp
            cp $hic_file inputs/hic/
        done

        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        mkdir verkko_hic

        all_wrapper.sh \
            assembly \
            verkko_hic \
            inputs

        tar -cvf ~{name}_~{tag}.tar verkko_hic
    >>>

    output {

        File hic_phased_run_folder = "~{name}_~{tag}.tar"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "humanpangenomics/verkko_hic@sha256:84ee7348dc954c9c28f78580f6551eaec19d85d813e656f275cde46133da0396"
        preemptible: preemptible
    }
}