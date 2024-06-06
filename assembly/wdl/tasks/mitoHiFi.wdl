version 1.0

import "../../../QC/wdl/tasks/extract_reads_toGZ.wdl" as extractReadsToGZ_t


workflow mitoHifi {

    input {
        File related_mito_fasta
        File related_mito_genbank
        Array[File] hifi_reads
        String sample_id

        ## optional filtering of hifi reads based on alignment to ref
        Boolean filter_hifi = true
        File? related_concat_mito_fasta

        ## for extracting reads to fastq.gz
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 512        
    }
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Calls [MitoHifi](https://github.com/marcelauliano/MitoHiFi/tree/master) to assemble representative mito contig."
    }

    parameter_meta {
        related_mito_fasta: "Reference mitochondria sequence (such as rCRS for humans) in fasta format."
        related_mito_genbank: "Reference mitochondria sequence (such as rCRS for humans) in genbank format."
        hifi_reads: "Hifi reads to generate mito contigs. Does not need to be filtered for mito sequences."
        filter_hifi: "Boolean for whether to filter reads based on alignment to reference. Default is true."
        related_concat_mito_fasta: "(Optional) Reference mitochondria sequence (such as rCRS for humans) concatenated to itself as one sequence in fasta format."
        referenceFasta: "(Optional) reference file for extracting hifi reads to fastq.gz"
    }

    scatter (readFile in hifi_reads) {
        call extractReadsToGZ_t.extractReadstoGZ as hifiExtractedGz {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB,
                dockerImage=dockerImage
        }

        if (filter_hifi) {
            call findMitoReads {
                input:
                    hifi_reads           = hifiExtractedGz.extractedRead,
                    concat_ref_mito      = related_concat_mito_fasta,
                    sample_id            = sample_id
            }
        }
    }


    ## use mito hifi reads if filtering was run
    hifi_for_assembly = select_first([findMitoReads.output, hifiExtractedGz.extractedRead])
    
    call mitoHifi {
        input:
            related_mito_fasta   = related_mito_fasta,
            related_mito_genbank = related_mito_genbank,
            hifi_reads           = hifi_for_assembly
    }

    output {
        File blastOutput       = blastFasta.blastOutput
        File parsedBlastOutput = parseBlastOutput.parsedBlastOutput
    }
}

task findMitoReads {

    input {
        File concat_ref_mito
        File hifi_reads
        String sample_id
        
        Int threadCount      = 32
        Int memSizeGB        = 64
        Int addldisk         = 64  
        String dockerImage   = "ghcr.io/marcelauliano/mitohifi@sha256:d8a72cac5add1d2836d9be4c2a50b32bb352699d1778d0ddcaaca79e1ce604ab" # master
    }

    Int read_size = ceil(size(hifi_reads, "GB"))
    Int final_disk_dize = 3*read_size + addldisk

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## align hifi reads to concatenated reference
        minimap2 \
            -t"~{threadCount}" \
            --secondary=no \
            -cx map-hifi \
            "~{concat_ref_mito}" \
            "~{sep('-r', hifi_reads)}" \
            -o "~{sample_id}_aligned_to_concat_mito.paf"

        ## select for reads that have alignment block length over threshold
        "~{sample_id}_aligned_to_concat_mito.paf"

        ## pull reads
        seqkit grep \
            -f reads_to_pull.txt \
            /private/groups/migalab/juklucas/mitohifi_test/HG02004/hifi/all_cat.fastq.gz \
            -o pulled_reads.fastq.gz
    >>>

    output {
        File blastOutput = blastOutputName
    }

    runtime {
        cpu: 
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task mitoHifi {

    input {
        File related_mito_fasta
        File related_mito_genbank
        File hifi_reads
        String sample_id 
        
        String kingdom       = "animal"
        Int perc_blast_match = 90
        Int genetic_code     = 2

        Int threadCount      = 24
        Int memSizeGB        = 64
        Int addldisk         = 64  
        String dockerImage   = "ghcr.io/marcelauliano/mitohifi@sha256:d8a72cac5add1d2836d9be4c2a50b32bb352699d1778d0ddcaaca79e1ce604ab" # master
    }
    
    parameter_meta {
        related_mito_fasta: "Reference mitochondria sequence (such as rCRS for humans) in fasta format."
        related_mito_genbank: "Reference mitochondria sequence (such as rCRS for humans) in genbank format."
        hifi_reads: "Hifi reads to generate mito contigs. Does not need to be filtered for mito sequences."
        sample_id: "ID to use for output naming"
        kingdom: "allowable options are : animal,plant,fungi"
        genetic_code: "Organism genetic code (following NCBI). Use 2 for the vertebrate mitochondrial code"
        perc_blast_match: "Percentage of query in the blast match with close-related mito. Recommend 90 for human."
    }

    Int read_size = ceil(size(hifi_reads), "GB"))
    Int final_disk_dize = 3*read_size + addldisk

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        mkdir ~{sample_id}
        cd ~{sample_id}

        mitohifi.py \
            "-r " ~{sep('-r', hifi_reads)} \
            -f ~{related_mito_fasta} \
            -g ~{related_mito_genbank} \
            -t ~{threadCount} \
            -a ~{kingdom} \
            -p ~{perc_blast_match} \
            -o ~{genetic_code}

        cd ..

        cp "~{sample_id}"/final_mitogenome.fasta "~{sample_id}_final_mitogenome.fasta"
        cp "~{sample_id}"/contigs_stats.tsv "~{sample_id}_contigs_stats.tsv"
        tar czvf "~{sample_id}_mito_hifi_working.tar.gz" ~{sample_id}/  

    >>>

    output {
        File mito_assembly       = "~{sample_id}_final_mitogenome.fasta"
        File mito_assembly_stats = "~{sample_id}_contigs_stats.tsv"
        File assembly_tar        = "~{sample_id}_mito_hifi_working.tar.gz"
    }

    runtime {
        cpu: 
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
