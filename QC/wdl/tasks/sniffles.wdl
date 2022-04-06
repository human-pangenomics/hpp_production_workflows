version 1.0

# This is a task level wdl workflow to run the program SNIFFLES

workflow runSniffles {
  input{
    File inputBam
    String outputName
    String dockerImage

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
      # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
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
