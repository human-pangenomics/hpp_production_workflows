version 1.0 

workflow ncbi_datasets_download_genome_wf {

    input {
        String genome_accession
        String sample_name
        String haplotype_string
        String output_tag        = "asm"
        Boolean reheader_fasta   = false    
    }

    call ncbi_datasets_download_genome {
        input:
            genome_accession  = genome_accession,
            sample_name       = sample_name,
            haplotype_string  = haplotype_string,
            output_tag        = output_tag,
            reheader_fasta    = reheader_fasta
    }
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Downloads genbank assembly using genbank genome accession. gzips and renames file."
    }

    parameter_meta {
        genome_accession: "Genbank (or RefSeq) genome accession. Example: GCA_042027545.1"
        sample_name: "Sample ID to use in output file name. Example: HG00408. Output will be named {sample_name}_{haplotype_string}_{output_tag}.fa.gz"
        haplotype_string: "Haplotype id to use in output file name. Example: mat/pat/hap1/hap2"
        output_tag: "Assembly info to put in fasta file name"
        reheader_fasta: "(default false) rename sequences with #sample_name#haplotype_string before sequence ID"
    }
    output {
        File genbank_fasta_gz     = ncbi_datasets_download_genome.fasta_gz
        File genbank_fasta_gz_gzi = ncbi_datasets_download_genome.fasta_gz_gzi
        File genbank_fasta_gz_fai = ncbi_datasets_download_genome.fasta_gz_fai
    }

}
task ncbi_datasets_download_genome {
    input {
        String genome_accession
        String sample_name
        String haplotype_string
        String output_tag
        Boolean reheader_fasta   

        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=64
        String dockerImage="humanpangenomics/ncbi-datasets-cli@sha256:c3b3bb734d2b7eb9cee874bb25ae94d57b8e11f8e06a0e090b7010d9193d37e0" # 16.31.0
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
            sed "s/^>\([^ ]*\).*/>~{sample_name}\#~{haplotype_string}\#\1/" ncbi_dataset/data/~{genome_accession}/*.fna \
                > ~{output_fasta_prefix}
        else
            ## copy and rename file
            cp ncbi_dataset/data/~{genome_accession}/*.fna ~{output_fasta_prefix}.fa
        fi

        ## compress and create indices for fasta
        bgzip -i ~{output_fasta_prefix}.fa
        samtools faidx ~{output_fasta_prefix}.fa.gz

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
    }
}

