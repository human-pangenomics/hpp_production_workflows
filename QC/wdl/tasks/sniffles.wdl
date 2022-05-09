version 1.0

# This is a task level wdl workflow to run the program SNIFFLES

workflow runSniffles {

  call Sniffles
  
  output{
    File outputFile = Sniffles.outputFile
  }
}

task Sniffles{
  input{
    File inputBam
    String outputName
    
    String dockerImage = "quay.io/biocontainers/sniffles:1.0.12--h8b12597_1"
    Int memSizeGB = 64
    Int threadCount = 32
    Int diskSizeGB = 64

  }
  command <<<
      # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        sniffles -m ~{inputBam} -v ~{outputName} -d 500 -n -1 -s 3

  >>>
  output{
    File outputFile = glob("*.vcf")[0]
  }
  runtime{
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSizeGB + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
