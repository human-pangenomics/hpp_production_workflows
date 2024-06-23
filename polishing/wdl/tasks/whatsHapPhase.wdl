version 1.0

# This is a task level wdl workflow to apply a set of variants to an assembly for polishing using bcftools consensus

workflow runWhatsHapPhase {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Phase variants in a vcf using a bamfile"
    }
    call WhatsHapPhase
    output {
        File phasedVcf = WhatsHapPhase.phasedVcf
    }
}

task WhatsHapPhase {
    input {
        File vcfFile
        File vcfFileIdx
        File refFile
        File refFileIdx
        File bamFile
        File bamFileIdx
        String outPrefix

        String dockerImage = "tpesout/whatshap:latest"
        Int memSizeGB = 128
        Int threadCount =16
        Int diskSizeGB = 256
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        # soft link data and indexes so they are in same place
        REF=$(basename ~{refFile})
        REF_IDX=$(basename ~{refFileIdx})

        ln -s ~{refFile} ./$REF
        ln -s ~{refFileIdx} ./$REF_IDX

        VCF=$(basename ~{vcfFile})
        VCF_IDX=$(basename ~{vcfFileIdx})

        ln -s ~{vcfFile} ./$VCF
        ln -s ~{vcfFileIdx} ./$VCF_IDX

        BAM=$(basename ~{bamFile})
        BAM_IDX=$(basename ~{bamFileIdx})

        ln -s ~{bamFile} ./$BAM
        ln -s ~{bamFileIdx} ./$BAM_IDX

        whatshap phase -o ~{outPrefix}.vcf.gz --indels -r $REF $VCF $BAM --ignore-read-groups
    >>>
    output {
        File phasedVcf = "~{outPrefix}.vcf.gz"
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
