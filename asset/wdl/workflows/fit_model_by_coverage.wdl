version 1.0

import "../tasks/cov2counts.wdl" as cov2counts_t
import "../tasks/fit_model.wdl" as fit_model_t

workflow fitModelBamInputOntHiFi{
    input {
        String sampleName
        String sampleSuffix
        File covFilesTarGz
    }

    call extractCovFiles {
        input:
            sampleName = sampleName,
            sampleSuffix = sampleSuffix,
            covFilesTarGz = covFilesTarGz
    }

    call cov2counts_t.cov2counts as cov2countsOnt{
        input:
            sampleName = sampleName,
            sampleSuffix = "${sampleSuffix}.ont",
            coverageFile = extractCovFiles.ontCoverage
    }

    call cov2counts_t.cov2counts as cov2countsHiFi{
        input:
            sampleName = sampleName,
            sampleSuffix = "${sampleSuffix}.hifi",
            coverageFile = extractCovFiles.hifiCoverage
    }

    call fit_model_t.fitModel as fitModelOnt{
        input:
            counts = cov2countsOnt.counts,
            sampleName = "${sampleName}",
            sampleSuffix = "${sampleSuffix}.ont"
    }

    call fit_model_t.fitModel as fitModelHiFi{
        input:
            counts = cov2countsHiFi.counts,
            sampleName = "${sampleName}",
            sampleSuffix = "${sampleSuffix}.hifi"
    }

    output{
        File hifiCoverageGz = cov2countsHiFi.coverageGz
        File ontCoverageGz = cov2countsOnt.coverageGz
        File hifiCounts = cov2countsHiFi.counts
        File ontCounts = cov2countsOnt.counts
        File hifiProbTable = fitModelHiFi.probabilityTable
        File ontProbTable = fitModelOnt.probabilityTable
    }
}

task extractCovFiles{
    input {
        String sampleName
        String sampleSuffix
        File covFilesTarGz
        # runtime configurations
        Int memSizeGB=4
        Int threadCount=2
        Int diskSizeGB=64
        Int preemptible=2
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
    }
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

        tar -xzf ~{covFilesTarGz} --strip-components 1
    >>>
    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
    }

    output{
        File ontCoverage = glob("*.ont.cov")[0]
        File hifiCoverage = glob("*.hifi.cov")[0]
    }
}

