version 1.0

# This is a task level wdl workflow to run a python script to filter SV calls

workflow runFilterSV {
    input {
        File inputVcf
        String outputName
        String? dockerImage

    }
    call Filter {
        input:
            inputVcf = inputVcf,
            outputName = outputName,
            dockerImage = dockerImage
    }
    output {
        File outputFile = Filter.outputFile
    }
}

task Filter{
    input {
        File inputVcf
        String outputName
        String dockerImage = "kishwars/t2t_polishing:0.1"
        Int memSizeGB = 64
        Int threadCount = 32
        Int diskSizeGB = 64
    }
    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ln -s /opt/filter.py .
        python3 filter.py ~{inputVcf} > ~{outputName}

    >>>
    output {
        File outputFile = glob("*filtered.vcf")[0]
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
