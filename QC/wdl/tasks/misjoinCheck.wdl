version 1.0

import "minigraph.wdl" as minigraph_wf

workflow misjoinCheck {
    input {
        File asm_fasta
        File ref_fasta
        File ref_centromere_bed
        String name
    }

    call minigraph_wf.minigraphMap as minigraph_t {
        input:
            inputFasta     = asm_fasta,
            inputGraph     = ref_fasta,
            sampleName     = name,
            args           = "-xasm"
    }

    call paftoolsMisjoinCheck {
        input:
            minigraphPaf   = minigraph_t.outputPafGZ,
            centromereBed  = ref_centromere_bed,
            sampleName     = name
    }

    output {
        File misjoinSummary = paftoolsMisjoinCheck.misjoinSummary
    }
}


task paftoolsMisjoinCheck {

    input {
        File minigraphPaf
        File centromereBed
        String sampleName
        String outputFileTag = "misjoin_summary"
        String args = "-ec"
        
        Int memSizeGB   = 4
        Int diskSizeGB  = 32
        String dockerImage = "humanpangenomics/hpp_paftools@sha256:3400cc89dafad584a5678735b1a6f039e3e8bb7e0f91bcc7610b1ab2ebea08a7" # minimap2 v2.28
    }

    String misjoinSummaryFn = "${sampleName}.${outputFileTag}.txt"

    command <<<

        set -eux -o pipefail
        
        ## paftools call to identify misjoins
        paftools.js misjoin ~{args} ~{centromereBed} ~{minigraphPaf} > ~{misjoinSummaryFn}
        
    >>>

    output {

        File misjoinSummary  = misjoinSummaryFn
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}