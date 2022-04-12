version 1.0

# This is a task level wdl workflow to run IRIS to filter HIFI SV calls

workflow runIris {
    input{
        File genomeIn
        File readsIn
        File vcfIn
        String vcfOut
        String IrisOut
        String? dockerImage

    }
    call Iris{
        input:
            genomeIn = genomeIn,
            readsIn = readsIn,
            vcfIn = vcfIn,
            vcfOut = vcfOut,
            IrisOut = IrisOut,
            dockerImage = dockerImage
    }
    output{
        File outputFile = Iris.outputFile
    }
}

task Iris{
    input{
        File genomeIn
        File readsIn
        File vcfIn
        String vcfOut
        String IrisOut
        String dockerImage = "quay.io/biocontainers/irissv:1.0.4--hdfd78af_2"
        Int memSizeGB = 64
        Int threadCount = 32
        Int diskSizeGB = 64

    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        iris --hifi --keep_long_variants --keep_files genome_in=~{genomeIn} reads_in=~{readsIn} vcf_in=~{vcfIn} vcf_out=~{vcfOut} out_dir=~{IrisOut}

    >>>
    output{
        File outputFile = glob("*.iris.vcf")[0]
    }
    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
