version 1.0

# This is a task level wdl to run assembly to assembly alignments with minimap2

workflow asm2asm_aligner_wf {

    call asm2asm_aligner
  
    output {
        File asm2asm_paf = asm2asm_aligner.paf
    }

    meta {
        description: "assembly to assembly alignment with minimap2. Outputs paf"
    }
}

task asm2asm_aligner {
    input{
        File asm_fasta
        File ref_fasta
        String preset      = "asm5"
        String otherArgs   = ""

        Int threadCount    = 32   
        Int memSizeGB      = 64
        Int diskSizeGB     = 64    
        String dockerImage = "quay.io/biocontainers/minimap2@sha256:5c8733f1483e5e466e6e40f155aa89591b79f64aac42e4bab28859f9afe22143" #2.27--he4a0461_1
    }

    parameter_meta {
        asm_fasta: "Query assembly to align against reference. Should be fa, fasta, fa, or fa.gz."
        ref_fasta: "Target/reference assembly. Could be anything, if you think about it. Should be fa, fasta, fa, or fa.gz, though."
        preset: "(default is asm5) Use asm5 for intra-species alignments with less than 1 percent divergence."
        otherArgs: "(default is empty string) Arguments to be passed to minimap2. For example '--cs' if the output will be passed to paftools.js"
    }
  
    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        ASM_NAME=$(basename "~{asm_fasta}" | sed -E 's/\.fa(\.gz)?$//;s/\.fasta(\.gz)?$//')
        REF_NAME=$(basename "~{ref_fasta}" | sed -E 's/\.fa(\.gz)?$//;s/\.fasta(\.gz)?$//')

        PAF_NAME="${ASM_NAME}_on_${REF_NAME}.paf"

        minimap2 \
            -cx \
            ~{preset} \
            -t~{threadCount} \
            ~{ref_fasta} \
            ~{asm_fasta} \
            ~{otherArgs} \
            > "$PAF_NAME"
  >>>

  output {
    File paf = glob("*.paf")[0]
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSizeGB + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}