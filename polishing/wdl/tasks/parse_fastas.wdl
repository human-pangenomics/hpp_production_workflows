  version 1.0

workflow runParseFastas {
    call parseFastas
    output {
      File Hap1RawFasta = parseFastas.hap1RawFasta
      File Hap1RawFastaIndex = parseFastas.hap1RawFastaIndex
      File Hap2RawFasta = parseFastas.hap2RawFasta
      File Hap2RawFastaIndex = parseFastas.hap2RawFastaIndex
      File dipRawFastaGz = parseFastas.dipRawFastaGz
      File dipRawFastaGzIndex = parseFastas.dipRawFastaGzIndex
    }
}

task parseFastas {
    input {
        File hap1Fasta
        File hap2Fasta
        String sampleName
        Int memSizeGB = 4
        Int threadCount = 1
        Int diskSizeGB = 48
        String dockerImage = "kishwars/pepper_deepvariant:r0.8"
    }


    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        HAP1_FILENAME=$(basename -- "~{hap1Fasta}")
        HAP2_FILENAME=$(basename -- "~{hap2Fasta}")

        mkdir hap1_output
        mkdir hap2_output

        if [[ $HAP1_FILENAME =~ \.gz$ ]]; then
            cp ~{hap1Fasta} hap1_output/~{sampleName}.hap1.fasta.gz
            gunzip -f hap1_output/~{sampleName}.hap1.fasta.gz
        else
            cp ~{hap1Fasta} hap1_output/~{sampleName}.hap1.fasta
        fi

        samtools faidx hap1_output/~{sampleName}.hap1.fasta

        if [[ $HAP2_FILENAME =~ \.gz$ ]]; then
            cp ~{hap2Fasta} hap2_output/~{sampleName}.hap2.fasta.gz
            gunzip -f hap2_output/~{sampleName}.hap2.fasta.gz
        else
            cp ~{hap2Fasta} hap2_output/~{sampleName}.hap2.fasta
        fi

        samtools faidx hap2_output/~{sampleName}.hap2.fasta

        cat hap1_output/~{sampleName}.hap1.fasta \
            hap2_output/~{sampleName}.hap2.fasta  \
            > ~{sampleName}.diploid.fasta

        bgzip ~{sampleName}.diploid.fasta

        samtools faidx ~{sampleName}.diploid.fasta.gz

    >>>

    output {
        File hap1RawFasta = glob("hap1_output/*hap1.fasta")[0]
        File hap1RawFastaIndex = glob("hap1_output/*hap1.fasta.fai")[0]
        File hap2RawFasta = glob("hap2_output/*hap2.fasta")[0]
        File hap2RawFastaIndex = glob("hap2_output/*hap2.fasta.fai")[0]
        File dipRawFastaGz = glob("*.diploid.fasta.gz")[0]
        File dipRawFastaGzIndex = glob("*.diploid.fasta.gz.fai")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
