version 1.0

workflow assign_chromosomes {
    input {
        File ref_asm
        File input_asm
        String assembly_name
        File chrom_assignment_ignore_bed
    }

    call mashmap_for_chrom_assignment {
        input:
            ref_asm       = ref_asm,
            input_asm     = input_asm,
            assembly_name = assembly_name
    }

    call find_gaps_and_telomeres {
        input:
            input_asm     = input_asm,
            assembly_name = assembly_name
    }

    call assign_chrom {
        input:
            mashmap_paf   = mashmap_for_chrom_assignment.mashmap_paf,
            telo_bed      = find_gaps_and_telomeres.telo_bed,
            gaps_bed      = find_gaps_and_telomeres.gaps_bed,
            ignore_bed    = chrom_assignment_ignore_bed,
            assembly_name = assembly_name
    }

    output {
        File mashmap_paf     = mashmap_for_chrom_assignment.mashmap_paf
        File telo_bed        = find_gaps_and_telomeres.telo_bed
        File gaps_bed        = find_gaps_and_telomeres.gaps_bed
        File chrom_alias_txt = assign_chrom.chrom_alias_txt
        File t2t_chrom_tsv   = assign_chrom.t2t_chrom_tsv
    }

    parameter_meta {
        ref_asm: "Reference assembly FASTA file. Recommended CHM13v2.0. Can be gzipped."
        input_asm: "Input assembly FASTA file. Can be gzipped."
        assembly_name: "Name of the assembly for output files"
        chrom_assignment_ignore_bed: "Regions to ignore during chromosome assignment. For example acro P arms and PARs."
    }

    meta {
        author: "Julian Lucas"
        description: "Workflow for chromosome assignment using mashmap. Outputs chromAlias and T2T assignment."
    }
}

task mashmap_for_chrom_assignment {
    input {
        File ref_asm
        File input_asm
        String assembly_name

        Int perc_ident      = 95

        Int thread_count    = 8
        Int mem_gb          = 32
        Int disk_size       = 64
    }

    command <<<
        set -eux -o pipefail

        # Handle potential gzipped inputs
        if [[ ~{ref_asm} == *.gz ]]; then
            gunzip -c ~{ref_asm} > ref.fa
            REF_ASM="ref.fa"
        else
            REF_ASM="~{ref_asm}"
        fi

        if [[ ~{input_asm} == *.gz ]]; then
            gunzip -c ~{input_asm} > input.fa
            INPUT_ASM="input.fa"
        else
            INPUT_ASM="~{input_asm}"
        fi

        # Run mashmap
        mashmap \
            --threads ~{thread_count} \
            --perc_identity ~{perc_ident} \
            --noSplit \
            -r ${REF_ASM} \
            -q ${INPUT_ASM} \
            -o mashmap_~{assembly_name}_pi~{perc_ident}.paf
    >>>

    output {
        File mashmap_paf = "mashmap_~{assembly_name}_pi~{perc_ident}.paf"
    }

    runtime {
        docker: "quay.io/biocontainers/mashmap@sha256:583811458b1fdee13028e5dc81573aa7edc0b18d285fa0cb2029496632ab5876"
        cpu: thread_count
        memory: "~{mem_gb}GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
    }
}

task find_gaps_and_telomeres {
    input {
        File input_asm
        String assembly_name

        Int thread_count    = 2
        Int mem_gb          = 16
        Int disk_size       = 64
    }

    command <<<
        set -eux -o pipefail

        # Handle potential gzipped input
        if [[ ~{input_asm} == *.gz ]]; then
            gunzip -c ~{input_asm} > input.fa
            INPUT_ASM="input.fa"
        else
            INPUT_ASM="~{input_asm}"
        fi

        # Find telomeres
        seqtk telo -s 10 ${INPUT_ASM} > ~{assembly_name}.telo.bed

        # Find gaps
        seqtk gap -l 2 ${INPUT_ASM} > ~{assembly_name}.gaps.bed
    >>>

    output {
        File telo_bed = "~{assembly_name}.telo.bed"
        File gaps_bed = "~{assembly_name}.gaps.bed"
    }

    runtime {
        docker: "quay.io/biocontainers/seqtk:1.4--he4a0461_2"
        cpu: thread_count
        memory: "~{mem_gb}GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
    }
}

task assign_chrom {
    input {
        File mashmap_paf
        File telo_bed
        File gaps_bed
        File ignore_bed
        String assembly_name

        Int thread_count    = 2
        Int mem_gb          = 16
        Int disk_size       = 64
    }

    command <<<
        set -eux -o pipefail

        chrom_alias_and_assignment.py \
            --mashmap ~{mashmap_paf} \
            --telomeres ~{telo_bed} \
            --gaps ~{gaps_bed} \
            --ignore ~{ignore_bed} \
            --out-prefix ~{assembly_name}
    >>>

    output {
        File chrom_alias_txt = "~{assembly_name}.chromAlias.txt"
        File t2t_chrom_tsv = "~{assembly_name}.t2t_chromosomes.tsv"
    }

    runtime {
        docker: "humanpangenomics/chrom_assignment@sha256:21c3f085320fe246d1fa862236a74440ffe02d09d42935408b48645c559c161d"
        cpu: thread_count
        memory: "~{mem_gb}GB"
        disks: "local-disk ~{disk_size} SSD"
        preemptible: 1
    }
}