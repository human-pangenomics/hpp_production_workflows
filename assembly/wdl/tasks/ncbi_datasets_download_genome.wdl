version 1.0 

workflow ncbi_datasets_download_genome_wf {

    input {
        String genome_accession
        String sample_name
        String haplotype_string    
    }

    call ncbi_datasets_download_genome {
        input:
            genome_accession  = genome_accession,
            sample_name       = sample_name,
            haplotype_string  = haplotype_string
    }
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Downloads genbank assembly using genbank genome accession. gzips and renames file."
    }

    parameter_meta {
        genome_accession: "Genbank (or RefSeq) genome accession. Example: GCA_042027545.1"
        sample_name: "Sample ID to use in output file name. Example: HG00408. Output will be named {sample_name}_{haplotype_string}_fa.gz"
        haplotype_string: "Haplotype id to use in output file name. Example: mat/pat/hap1/hap2"
    }
    output {
        File genbank_fa_gz = ncbi_datasets_download_genome.fasta_gz
    }

}
task ncbi_datasets_download_genome {
    input {
        String genome_accession
        String sample_name
        String haplotype_string

        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=64
        String dockerImage="biocontainers/ncbi-datasets-cli@sha256:3535d6bc909dbee0e72abfa9ca079a422b08b03dced4bbf7e0affae2a0118e9e" # 16.22.1_cv1 (amd64)
        Int preemptible=2
    }

    String output_fasta_name = "~{sample_name}_~{haplotype_string}_fa.gz"

    command <<<
        # Stop on any command failure, return failure if any command in a pipeline fails.
        # Prevent using unset variables.
        set -euo pipefail
        ## print commands before running
        set -o xtrace
        
        datasets download genome accession ~{genome_accession} --filename "~{genome_accession}_dataset.zip"

        ## unzip archive
        unzip "~{genome_accession}_dataset.zip"

        gzip -c "ncbi_dataset/data/~{genome_accession}/*.fna" > "~{output_fasta_name}"

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fasta_gz = output_fasta_name
    }
}

