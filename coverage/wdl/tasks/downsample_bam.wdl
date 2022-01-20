version 1.0

workflow runDownsampleBam{
    call downsampleBam
    output{
        File downsampledBam = downsampleBam.downsampledBam
        File downsampledBai = downsampleBam.downsampledBai
    }
}

task downsampleBam{
    input{
        File bam
        Float rate
        # runtime configurations
        Int memSize=32
        Int threadCount=8
        Int diskSize=512
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
        String zones
    }
    command <<<
        set -eu -o pipefail -o xtrace
        
        FILENAME=$(basename ~{bam})
        PREFIX=${FILENAME%%.bam}

        mkdir output
        java -jar /home/apps/picard.jar DownsampleSam \
                                        I=~{bam} \
                                        O=output/${PREFIX}.downsampled.bam \
                                        P=~{rate}
        samtools index -@ ~{threadCount} output/${PREFIX}.downsampled.bam
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
    output {
        File downsampledBam = glob("output/*.bam")[0]
        File downsampledBai = glob("output/*.bai")[0]
    }
}

