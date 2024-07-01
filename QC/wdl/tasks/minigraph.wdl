version 1.0

workflow runMinigraphMap {

    call minigraphMap 

    output {
        File outputPafGZ = minigraphMap.outputPafGZ
    }
}



task minigraphMap {

    input {
        File inputFasta
        File inputGraph
        String sampleName
        String outputFileTag = "minigraph"
        String args = ""
        
        Int threadCount = 8
        Int memSizeGB   = 64
        Int diskSizeGB  = 64
        String dockerImage = "humanpangenomics/hpp_minigraph@sha256:51cbfe78ca8b9c47de415c2f72a534c6124cd8224302fb2ae3d6e4707c2618df" # v0.21
    }

    String outputPafGz = "${sampleName}.${outputFileTag}.paf.gz"

    command <<<

        set -eux -o pipefail
        
        ## Run minigraph
        minigraph \
            -t ~{threadCount} \
            ~{args} \
            ~{inputGraph} \
            ~{inputFasta} \
            | gzip > ~{outputPafGz}

    >>>

    output {

        File outputPafGZ  = outputPafGz
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}