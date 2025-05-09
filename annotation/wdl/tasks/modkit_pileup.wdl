version 1.0

workflow hprc_modkit_pileup {

    meta {
        author: "Julian Menendez"
        email: "jmmenend@ucsc.edu"
        description: "Runs modkit pileup to produce bigwig of methylation. Requires assembly and diploid aligned bam with MM/ML tags."
    }

    parameter_meta {
        sample_id: "Sample identifier used for naming intermediate files. Not included in output files."
        hap1_prefix: "Prefix to use for naming haplotype 1 output files (example 'HG002_hap1_assembly')"
        hap2_prefix: "Prefix to use for naming haplotype 2 output files (example 'HG002_hap2_assembly')"
        output_suffix: "String suffix to append to output filenames for identification (example: 'modkit_ont')"
        hap1_fasta: "FASTA file containing the first haplotype assembly (allowed formats: .fasta, .fa, .fasta.gz, .fa.gz)"
        hap2_fasta: "FASTA file containing the second haplotype assembly (allowed formats: .fasta, .fa, .fasta.gz, .fa.gz)"
        methyl_bam: "BAM file of reads with methylation tags (MM/ML) aligned to diploid assembly"
        extra_samtools_flags: "Flags to pass to samtools during BAM filtering (default: '-q 10 -m 10000' for min mapq 10, min aligned length 10kb)"
        modkit_thresholds: "Thresholds for modkit pileup to filter methylation calls (default for ONT: '--mod-thresholds m:0.8 --filter-threshold C:0.5')"
    }

    input {
        String sample_id

        String output_suffix

        File hap1_fasta
        File hap2_fasta
        String hap1_prefix
        String hap2_prefix

        File methyl_bam

        String extra_samtools_flags = "-q 10 -m 10000" # min mapq 10, min alen 10kb
        String modkit_thresholds = "--mod-thresholds m:0.8 --filter-threshold C:0.5"
    }

    Boolean hap2_is_gzipped = sub(basename(hap2_fasta), "^.+\\.", "") == "gz"
    if (hap2_is_gzipped) {
        call unzip_fasta as unzip_hap2 {
            input:
                ref_fasta_gz = hap2_fasta
        }
    }

    Boolean hap1_is_gzipped = sub(basename(hap1_fasta), "^.+\\.", "") == "gz"
    if (hap1_is_gzipped) {
        call unzip_fasta as unzip_hap1 {
            input:
                ref_fasta_gz = hap1_fasta
        }
    }

    call index_fa as index_hap1 {
        input:
            fasta = select_first([unzip_hap1.fasta, hap1_fasta])
    }

    call index_fa as index_hap2 {
        input:
            fasta = select_first([unzip_hap2.fasta, hap2_fasta])
    }

    call merge_haps {
        input:
            hap1_fasta = index_hap1.fa,
            hap2_fasta = index_hap2.fa,
            sample_id = sample_id
    }

    call index_fa as index_diploid {
        input:
            fasta = merge_haps.diploid_fa
    }

    call index_bam {
        input:
            input_bam = methyl_bam,
            extra_samtools_flags = extra_samtools_flags,
            sample_id = sample_id,
            output_suffix = output_suffix
    }

    call modkit_pileup {
        input:
            bam = index_bam.bam,
            bai = index_bam.bam_bai,
            ref_fa = index_diploid.fa,
            ref_fa_fai = index_diploid.fa_fai,
            modkit_thresholds = modkit_thresholds,
            output_suffix = output_suffix,
            sample_id = sample_id
    }

    call split_beds_by_hap {
        input:
            pileup = modkit_pileup.modkit_pileup_bed,
            pileup_bedgraph = modkit_pileup.modkit_pileup_bedgraph,
            hap1_chrom_names = merge_haps.hap1_chrom_names,
            hap1_prefix = hap1_prefix,
            hap2_chrom_names = merge_haps.hap2_chrom_names,
            hap2_prefix = hap2_prefix,
            output_suffix = output_suffix
    }

    call bedgraph_to_bigwig as hap1_bigwig {
        input:
            pileup_bedgraph = split_beds_by_hap.hap1_pileup_bedgraph,
            fa_genome = index_hap1.fa_genome
    }

    call bedgraph_to_bigwig as hap2_bigwig {
        input:
            pileup_bedgraph = split_beds_by_hap.hap2_pileup_bedgraph,
            fa_genome = index_hap2.fa_genome
    }

    output {
        File hap1_pileup_bed = split_beds_by_hap.hap1_pileup
        File hap1_pileup_bigwig = hap1_bigwig.pileup_bigwig

        File hap2_pileup_bed = split_beds_by_hap.hap2_pileup
        File hap2_pileup_bigwig = hap2_bigwig.pileup_bigwig
    }
}

task unzip_fasta {
    input {
        File ref_fasta_gz

        Int threadCount = 4
        Int memSizeGB  = 16
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(ref_fasta_gz, "GB"))       
    Int final_disk_size = input_file_size * 9

    String fa_basename = basename(ref_fasta_gz)
    String fasta_output = sub(fa_basename, "\\.gz$", "")

    command <<<
        set -eux -o pipefail

        gzip -c -d ~{ref_fasta_gz} > ~{fasta_output}
    >>>

    output {
        File fasta = fasta_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "ubuntu:24.04"
        preemptible: 1
    }
}

task index_fa {
    input {
        File fasta

        Int threadCount = 4
        Int memSizeGB = 16
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(fasta, "GB"))       
    Int final_disk_size = input_file_size * 3

    String fa_basename = basename(fasta)
    String fa_index_output = "~{fa_basename}.fai"
    String fa_genome_output = "~{fa_basename}.genome"

    command <<<
        set -eux -o pipefail

        samtools faidx --fai-idx ~{fa_index_output} ~{fasta}
        cut -f 1,2 ~{fa_index_output} > ~{fa_genome_output}
    >>>

    output {
        File fa = fasta
        File fa_fai = fa_index_output
        File fa_genome = fa_genome_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "quay.io/staphb/samtools:1.21"
        preemptible: 1
    }
}

task merge_haps {
    input {
        File hap1_fasta
        File hap2_fasta

        String sample_id

        Int threadCount = 8
        Int memSizeGB = 32
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(hap1_fasta, "GB")) + ceil(size(hap2_fasta, "GB"))       
    Int final_disk_size = input_file_size * 6

    String merged_fasta_output = "~{sample_id}_merged.fasta"

    String hap1_chrom_names_output = "~{sample_id}.hap1_chrom_names.txt"
    String hap2_chrom_names_output = "~{sample_id}.hap2_chrom_names.txt"

    command <<<
        set -eux -o pipefail

        cat ~{hap2_fasta} ~{hap1_fasta} > ~{merged_fasta_output}

        # extract the chromosome names from hap2_fasta and hap1_fasta and make them outputs
        grep "^>" ~{hap1_fasta} | sed 's/^>//' > ~{hap1_chrom_names_output}
        grep "^>" ~{hap2_fasta} | sed 's/^>//' > ~{hap2_chrom_names_output}
    >>>

    output {
        File diploid_fa = merged_fasta_output

        File hap1_chrom_names = hap1_chrom_names_output
        File hap2_chrom_names = hap2_chrom_names_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "ubuntu:24.04"
        preemptible: 1
    }
}

task index_bam {
    input {
        File input_bam
        String sample_id

        String output_suffix

        String extra_samtools_flags

        Int threadCount = 16
        Int memSizeGB = 128
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(input_bam, "GB"))       
    Int final_disk_size = input_file_size * 3

    String bam_output = "~{sample_id}.~{output_suffix}.bam"
    String index_output = "~{sample_id}.~{output_suffix}.bam.bai"

    command <<<
        set -eux -o pipefail

        # samtools view/filter and sort
        samtools view -bh ~{extra_samtools_flags} -@ ~{threadCount} ~{input_bam} | \
            samtools sort -@ ~{threadCount} -o ~{bam_output}

        # samtools index
        samtools index -@ ~{threadCount} -o ~{index_output} ~{bam_output}
    >>>

    output {
        File bam = bam_output
        File bam_bai = index_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "quay.io/staphb/samtools:1.21"
        preemptible: 1
    }
}

task modkit_pileup {
    input {
        File ref_fa 
        File ref_fa_fai

        File bam
        File bai

        String sample_id
        String output_suffix

        String modkit_thresholds

        Int threadCount = 16
        Int memSizeGB = 128
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(bam, "GB"))       
    Int final_disk_size = input_file_size * 3

    String modkit_pileup_bed_output = "~{sample_id}.~{output_suffix}.CpG.pileup.bed"
    String modkit_pileup_bedgraph_output = "~{sample_id}.~{output_suffix}.5mC.bedgraph"

    String ref_basename = basename(ref_fa)
    String ref_fai_basename = basename(ref_fa_fai)

    command <<<
        set -eux -o pipefail

        ln -s ~{ref_fa} ./~{ref_basename}
        ln -s ~{ref_fa_fai} ./~{ref_fai_basename}

        modkit pileup \
            -t ~{threadCount} \
            --force-allow-implicit \
            --cpg \
            --ref ./~{ref_basename} \
            --combine-strands \
            ~{modkit_thresholds} \
            ~{bam} ~{modkit_pileup_bed_output}

        # use awk to create bedgraph
        awk -v OFS='\t' '{if ($4=="m") print $1, $2, $3, $11}' ~{modkit_pileup_bed_output} > ~{modkit_pileup_bedgraph_output}
    >>>

    output {
        File modkit_pileup_bed = modkit_pileup_bed_output
        File modkit_pileup_bedgraph = modkit_pileup_bedgraph_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "jmmenend/modkit:0.4.2"
        preemptible: 1
    }
}

task split_beds_by_hap {
    input {
        File pileup
        File pileup_bedgraph

        File hap1_chrom_names
        File hap2_chrom_names

        String hap1_prefix
        String hap2_prefix
        String output_suffix

        Int threadCount = 8
        Int memSizeGB  = 64
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(pileup, "GB")) + ceil(size(pileup_bedgraph, "GB"))       
    Int final_disk_size = input_file_size * 6

    String hap1_pileup_output = "~{hap1_prefix}.~{output_suffix}.CpG_pileup.bed"
    String hap1_pileup_bedgraph_output = "~{hap1_prefix}.~{output_suffix}.5mC.bedgraph"

    String hap2_pileup_output = "~{hap2_prefix}.~{output_suffix}.CpG_pileup.bed"
    String hap2_pileup_bedgraph_output = "~{hap2_prefix}.~{output_suffix}.5mC.bedgraph"

    command <<<
        set -eux -o pipefail

        # split pileups by haplotype
        grep -wf ~{hap1_chrom_names} ~{pileup} > ~{hap1_pileup_output}
        grep -wf ~{hap2_chrom_names} ~{pileup} > ~{hap2_pileup_output}

        # split bedgraphs by haplotype
        grep -wf ~{hap1_chrom_names} ~{pileup_bedgraph} > ~{hap1_pileup_bedgraph_output}
        grep -wf ~{hap2_chrom_names} ~{pileup_bedgraph} > ~{hap2_pileup_bedgraph_output}
    >>>

    output {
        File hap1_pileup = hap1_pileup_output
        File hap2_pileup = hap2_pileup_output

        File hap1_pileup_bedgraph = hap1_pileup_bedgraph_output
        File hap2_pileup_bedgraph = hap2_pileup_bedgraph_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "ubuntu:24.04"
        preemptible: 1
    }
}

task bedgraph_to_bigwig {
    input {
        File pileup_bedgraph
        File fa_genome

        Int threadCount = 16
        Int memSizeGB = 128
    }

    # Estimate disk size required
    Int input_file_size = ceil(size(pileup_bedgraph, "GB"))     
    Int final_disk_size = input_file_size * 6

    String bedgraph_basename = basename(pileup_bedgraph, ".bedgraph")
    String pileup_bigwig_output = "~{bedgraph_basename}.bigwig"

    command <<<
        set -eux -o pipefail

        bedGraphToBigWig \
            ~{pileup_bedgraph} \
            ~{fa_genome} \
            ~{pileup_bigwig_output}
    >>>

    output {
        File pileup_bigwig = pileup_bigwig_output
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_size + " SSD"
        docker: "jmmenend/bedgraphtobigwig:latest"
        preemptible: 1
    }
}