version 1.0

import "subBamByBed.wdl" as subBamByBed_t

workflow runGetMapQTable {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "make table of mapq values for correctBam"
    }
    input {
        File allHifiToHap2Bam
        File allHifiToHap2Bai
        File allHifiToHap1Bam
        File allHifiToHap1Bai
        File secPhaseBed
    }
    call subBamByBed_t.SubBamByBed as subBamHap1 {
        input:
            Bam=allHifiToHap1Bam,
            Bai=allHifiToHap1Bai,
            Bed=secPhaseBed
    }
    call subBamByBed_t.SubBamByBed as subBamHap2 {
        input:
            Bam=allHifiToHap2Bam,
            Bai=allHifiToHap2Bai,
            Bed=secPhaseBed
    }
    call getMapQTable {
        input:
            allHifiToHap2Bam=subBamHap2.subBam,
            allHifiToHap2Bai=subBamHap2.subBai,
            allHifiToHap1Bam=subBamHap1.subBam,
            allHifiToHap1Bai=subBamHap1.subBai
    }

    output {
        File mapqTable = getMapQTable.mapqTable
    }
}

task getMapQTable {
    input {
        File allHifiToHap2Bam
        File allHifiToHap2Bai
        File allHifiToHap1Bam
        File allHifiToHap1Bai

        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
        Int threads = 8
        Int memSizeGb = 128
        Int diskSizeGb = 256
    }
    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # softlink bam and index so they are in same directory

        ln -s ~{allHifiToHap2Bam} Hap2.bam
        ln -s ~{allHifiToHap2Bai} Hap2.bai

        ln -s ~{allHifiToHap1Bam} Hap1.bam
        ln -s ~{allHifiToHap1Bai} Hap1.bai

        # Set param file based on input hifi or ont read alignments
        samtools view Hap1.bam | awk '{print $1"\t"$3"\t"$4"\t"$5}' > mapq_table_Hap1.tsv
        samtools view Hap2.bam | awk '{print $1"\t"$3"\t"$4"\t"$5}' > mapq_table_Hap2.tsv
        cat mapq_table_Hap1.tsv mapq_table_Hap2.tsv > mapq_table_all.tsv

    >>>
    output {
        File mapqTable = "mapq_table_all.tsv"
    }

    runtime {
        memory: memSizeGb + " GB"
        cpu: threads
        disks: "local-disk " + diskSizeGb + " SSD"
        docker: dockerImage
    }
}
