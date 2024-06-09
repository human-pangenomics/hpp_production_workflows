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
        referenceFasta: "(Optional) reference file for extracting hifi reads to fastq.gz"
    }

    ## Extract reads to fastq.gz, filter by alignment block length if filter_hifi is true
    scatter (readFile in hifi_reads) {
        call extractReadsToGZ_t.extractReadstoGZ as hifiExtractedGz {
            input:
                readFile=readFile,
                referenceFasta=referenceFasta,
                memSizeGB=4,
                threadCount=4,
                diskSizeGB=fileExtractionDiskSizeGB
        }

        if (filter_hifi) {
            call findMitoReads {
                input:
                    hifi_fastq_gz  = hifiExtractedGz.extractedRead,
                    ref_mito       = related_mito_fasta,
                    sample_id      = sample_id
            }
        }
        File hifi_for_assembly = select_first([findMitoReads.filtered_reads, hifiExtractedGz.extractedRead])
    }
    
    ## assembly mito genome
    call mitoHifi {
        input:
            related_mito_fasta   = related_mito_fasta,
            related_mito_genbank = related_mito_genbank,
            hifi_reads           = hifi_for_assembly
    }

    ## align assembled genome against reference and call variants (to run in downstream QC)
    call mitoAssemblyVariants {
        input:
            related_mito_fasta   = related_mito_fasta,
            assembled_mito       = mitoHifi.mito_assembly,
            sample_id            = sample_id
    }


    output {
        ## MitoHiFi outputs
        File mitoHiFi_assembly     = mitoHifi.mito_assembly
        File mitoHiFi_stats        = mitoHifi.mito_assembly_stats
        File mitoHiFi_output_tar   = mitoHifi.assembly_tar

        ## mito assembly evaluations
        File mitoHiFi_eval_bam     = mitoAssemblyVariants.mito_asm_on_ref_bam
        File mitoHiFi_eval_vcf     = mitoAssemblyVariants.mito_asm_on_ref_vcf
    }
}

task findMitoReads {

    input {
        File ref_mito
        File hifi_fastq_gz
        String sample_id

        Int min_alignment_block_len = 12000
        
        Int threadCount      = 32
        Int memSizeGB        = 64
        Int addldisk         = 64  
        String dockerImage   = "humanpangenomics/mitohifi_utils@sha256:cbd31256f13772d442f80a2e37cd2f36703a780a852d59dd6028e2cf8892c7d7"
    }

    Int read_size = ceil(size(hifi_fastq_gz, "GB"))    
    Int final_disk_dize = 3*read_size + addldisk

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## concatenate reference to itself to help with problems caused by circular genomes
        input_mito="~{ref_mito}"

        python3 << EOF
        from Bio import SeqIO

        def self_concatenate_fasta_sequences(file1, output_file):
            seq1 = SeqIO.read(file1, "fasta")
            concatenated_sequence = seq1.seq + seq1.seq
            
            # Create a new SeqRecord with the concatenated sequence
            concatenated_record = SeqIO.SeqRecord(concatenated_sequence, id=f"{seq1.id}-self-concatenated", description=seq1.description)
            SeqIO.write(concatenated_record, output_file, "fasta")

        input_mito_py = "$input_mito"
        output_file   = "mito_self_concat.fasta"

        self_concatenate_fasta_sequences(input_mito_py, output_file)
        EOF


        ## align hifi reads to concatenated reference
        minimap2 \
            -t"~{threadCount}" \
            --secondary=no \
            -cx map-hifi \
            mito_self_concat.fasta \
            "~{hifi_fastq_gz}" \
            -o "~{sample_id}_aligned_to_concat_mito.paf"


        ## select for reads that have alignment block length over threshold
        min_len=~{min_alignment_block_len}

        awk -v min_len="$min_len" '$11 >= min_len {print $1}' \
            "~{sample_id}_aligned_to_concat_mito.paf" \
            > reads_to_pull.txt


        ## pull reads
        FILENAME=$(basename -- "~{hifi_fastq_gz}" | sed 's/.gz$//' ) 
        seqkit grep \
            -f reads_to_pull.txt \
            "~{hifi_fastq_gz}" \
            -o "${FILENAME}_filtered_over_~{min_alignment_block_len}.fastq.gz"
    >>>

    output {
        File read_aligned_to_concat_mito = "~{sample_id}_aligned_to_concat_mito.paf"
        File filtered_reads = glob("*_filtered_over_~{min_alignment_block_len}.fastq.gz")[0]
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task mitoHifi {
    
    parameter_meta {
        related_mito_fasta: "Reference mitochondria sequence (such as rCRS for humans) in fasta format."
        related_mito_genbank: "Reference mitochondria sequence (such as rCRS for humans) in genbank format."
        hifi_reads: "Hifi reads to generate mito contigs. Does not need to be filtered for mito sequences."
        sample_id: "ID to use for output naming"
        kingdom: "allowable options are : animal,plant,fungi"
        genetic_code: "Organism genetic code (following NCBI). Use 2 for the vertebrate mitochondrial code"
        perc_blast_match: "Percentage of query in the blast match with close-related mito. Recommend 90 for human."
    }

    input {
        File related_mito_fasta
        File related_mito_genbank
        Array[File] hifi_reads
        String sample_id 
        
        String kingdom       = "animal"
        Int perc_blast_match = 90
        Int genetic_code     = 2

        Int threadCount      = 24
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


        cat ~{sep=" " hifi_reads} > all_reads.fastq.gz

        mkdir ~{sample_id}
        cd ~{sample_id}

        mitohifi.py \
            -r  ../all_reads.fastq.gz \
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
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task mitoAssemblyVariants {

    input {
        File related_mito_fasta
        File assembled_mito
        String sample_id 

        Int threadCount      = 24
        Int memSizeGB        = 64
        Int diskSizeGB       = 64  
        String dockerImage   = "humanpangenomics/mitohifi_utils@sha256:cbd31256f13772d442f80a2e37cd2f36703a780a852d59dd6028e2cf8892c7d7"

    }

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        minimap2 \
            -t~{threadCount} \
            --secondary=no \
            -ax asm5 \
            ~{related_mito_fasta} \
            ~{assembled_mito} \
            | samtools view -Sb \
            | samtools sort \
            > "~{sample_id}_mito_asm_on_reference.bam"

        samtools index "~{sample_id}_mito_asm_on_reference.bam"

        bcftools mpileup \
            -f "~{related_mito_fasta}" \
            "~{sample_id}_mito_asm_on_reference.bam" \
            | bcftools call \
                -mv -Ov \
                --ploidy 1 \
                -o "~{sample_id}_mito_asm_on_reference.vcf"

    >>>

    output {
        File mito_asm_on_ref_bam = "~{sample_id}_mito_asm_on_reference.bam"
        File mito_asm_on_ref_vcf = "~{sample_id}_mito_asm_on_reference.vcf"
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
