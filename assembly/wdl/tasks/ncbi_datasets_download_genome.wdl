version 1.0 

workflow ncbi_datasets_download_genome_wf {

    input {
        String sample_name

        Int haplotype_int
        String haplotype_string
        String genome_accession         

        Int? haplotype_int_2
        String? haplotype_string_2
        String? genome_accession_2

        String output_tag        = "asm"
        Boolean reheader_fasta   = false    
    }

    ## download genbank assembly
    call ncbi_datasets_download_genome as download_asm1 {
        input:
            sample_name       = sample_name,
            genome_accession  = genome_accession,
            haplotype_int     = haplotype_int,
            haplotype_string  = haplotype_string,
            output_tag        = output_tag,
            reheader_fasta    = reheader_fasta
    }
    

    if (defined(genome_accession_2) && defined(haplotype_string_2)) {

        ## download second genbank assembly (other haplotype of the same sample)
        call ncbi_datasets_download_genome as download_asm2 {
            input:
                sample_name       = sample_name,
                genome_accession  = select_first([genome_accession_2]),
                haplotype_int     = select_first([haplotype_int_2]),
                haplotype_string  = select_first([haplotype_string_2]),
                output_tag        = output_tag,
                reheader_fasta    = reheader_fasta
        }
    }

    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Downloads genbank assembly using genbank genome accession. gzips and renames file. Optionally downloads two haplotypes."
    }

    parameter_meta {
        sample_name: "Sample ID to use in output file name. Example: HG00408. Output will be named {sample_name}_{haplotype_string}_{output_tag}.fa.gz"
        genome_accession: "Genbank (or RefSeq) genome accession. Example: GCA_042027545.1"        
        haplotype_int: "Haplotype integer to use in sequence names (according to PanSN spec). 1/2 is used for pat/mat, respectively."
        haplotype_string: "Haplotype id to use in output file name. Example: mat/pat/hap1/hap2"
        genome_accession_2: "Accession for second haplotype, if one is present."
        haplotype_string_2: "Haplotype id for second haplotype, if one is present."
        output_tag: "Assembly info to put in fasta file name"
        reheader_fasta: "(default false) rename sequences with #sample_name#haplotype_string before sequence ID"
    }

    output {
        File genbank_fa_hap1_gz      = download_asm1.fasta_gz
        File genbank_fa_hap1_gz_gzi  = download_asm1.fasta_gz_gzi
        File genbank_fa_hap1_gz_fai  = download_asm1.fasta_gz_fai
        File genbank_fa_hap1_gz_md5  = download_asm1.fasta_gz_md5
        
        File? genbank_fa_hap2_gz     = download_asm2.fasta_gz
        File? genbank_fa_hap2_gz_gzi = download_asm2.fasta_gz_gzi
        File? genbank_fa_hap2_gz_fai = download_asm2.fasta_gz_fai        
        File? genbank_fa_hap2_gz_md5 = download_asm2.fasta_gz_md5
    }

}
task ncbi_datasets_download_genome {
    input {
        String genome_accession
        String sample_name
        Int haplotype_int
        String haplotype_string
        String output_tag
        Boolean reheader_fasta   

        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=64
        String dockerImage="humanpangenomics/ncbi-datasets-cli@sha256:bad03e0c7e44f9e9038978d2a570b9d112d58d325f97c7bb6b9de4721b018a31" # 16.31.0
        Int preemptible=2
    }

    String output_fasta_prefix = "~{sample_name}_~{haplotype_string}_~{output_tag}"

    command <<<
        # Stop on any command failure, return failure if any command in a pipeline fails.
        # Prevent using unset variables.
        set -euo pipefail
        ## print commands before running
        set -o xtrace
        
        datasets download genome accession ~{genome_accession} --filename "~{genome_accession}_dataset.zip"

        ## unzip archive
        unzip "~{genome_accession}_dataset.zip"

        if [[ ~{reheader_fasta} == true ]]; then
            ## reheader file from: 
            ## CM086346.1 Homo sapiens isolate HG00609 chromosome 1, whole genome shotgun sequence
            ## to:
            ## HG00609#hap1#CM086352.1 (for example)
            sed "s/^>\([^ ]*\).*/>~{sample_name}\#~{haplotype_int}\#\1/" ncbi_dataset/data/~{genome_accession}/*.fna \
                > ~{output_fasta_prefix}.masked.fa
        else
            ## copy and rename file
            cp ncbi_dataset/data/~{genome_accession}/*.fna ~{output_fasta_prefix}.masked.fa
        fi

        ## genbank soft masks and has a line width of 80 characters.
        ## remove soft masking from Genbank and set line width to 60
        ## sort sequences
        seqkit seq -u -w 60 ~{output_fasta_prefix}.masked.fa \
            | seqkit sort --natural-order --two-pass \
            > ~{output_fasta_prefix}.fa


        ## compress and create indices for fasta
        bgzip -i ~{output_fasta_prefix}.fa
        samtools faidx ~{output_fasta_prefix}.fa.gz


        md5sum ~{output_fasta_prefix}.fa.gz > ~{output_fasta_prefix}.fa.gz.md5
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fasta_gz     = "~{output_fasta_prefix}.fa.gz"
        File fasta_gz_gzi = "~{output_fasta_prefix}.fa.gz.gzi"
        File fasta_gz_fai = "~{output_fasta_prefix}.fa.gz.fai"
        File fasta_gz_md5 = "~{output_fasta_prefix}.fa.gz.md5"
    }
}

