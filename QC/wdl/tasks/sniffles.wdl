version 1.0

# This is a task level wdl workflow to run the program SNIFFLES

workflow runSniffles {
  input{
    File inputBam
    String outputName
    String? dockerImage

  }
  call Sniffles{
    input:
      inputBam = inputBam,
      outputName = outputName,
      dockerImage = dockerImage
  }
  output{
    File outputFile = Sniffles.outputFile
  }
}

task Sniffles{
  input{
    File inputBam
    String outputName
    String dockerImage = "quay.io/biocontainers/sniffles:1.0.12--h8b12597_1"

  }
  # sniffles -m "${INPUT_DIR}"/chm13.draft_v0.9.hifi_20k.wm_1.11.pri.chr20.bam -v "${OUTPUT_DIR}"/"${HIFI_SVS_SNIFFLES}" -d 500 -n -1 -s 3
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
    docker: dockerImage
  }
}
