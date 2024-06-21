version 1.0

workflow runMarginPhase {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "phase variants with margin phase"
    }

    call marginPhase {
    }

    output {
        File out_margin_phase_svs = marginPhase.phasedVcf
    }
}

task marginPhase {
    input {
        File vcfFile
        File vcfFileIdx
        File refFile
        File refFileIdx
        File bamFile
        File bamFileIdx
        String outPrefix
        String HifiOrONT

        File? configFile # if config file is passed in, it will override HifiOrONT

        String dockerImage = "miramastoras/marginphase_sv:latest"
        Int threads = 32
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # Set param file based on input hifi or ont read alignments

        if [[ -f "~{configFile}" ]] ; then
            PARAMS=~{configFile}
        elif [[ ~{HifiOrONT} =~ Hifi ]]; then
            PARAMS=/opt/margin/params/phase/allParams.phase_vcf.pb-hifi.json
        else
            PARAMS=/opt/margin/params/phase/allParams.phase_vcf.ont.json
        fi

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

        mkdir output/
        margin phase ${BAM} ${REF} ${VCF} ${PARAMS} -t ~{threads} -o output/~{outPrefix} -M


    >>>
    output {
        File phasedVcf = glob("output/*.vcf")[0]
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
