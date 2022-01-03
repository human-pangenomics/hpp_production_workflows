version 1.0

import "../tasks/cov2counts.wdl" as cov2counts_t
import "../tasks/fit_model.wdl" as fit_model_t
import "../tasks/find_blocks.wdl" as find_blocks_t
import "../tasks/subset_coverage.wdl" as sub_cov_t

workflow runFitModelBed{
    input {
        File bed
        String suffix
        File coverageGz
        Float covFloat
    }
    if ("${bed}" == "") {
        call makeNullOutput{
            input:
                prefix = basename("${coverageGz}", ".cov.gz"),
                suffix = suffix
        }
    }
    if("${bed}" != ""){
        call sub_cov_t.subsetCoverage{
            input:
                coverageGz = coverageGz,
                blocksBed = bed,
                suffix = suffix
        }
        call cov2counts_t.cov2counts{
            input:
                coverageGz = subsetCoverage.outputCoverageGz
        }
        call fit_model_t.fitModel {
            input:
                counts = cov2counts.counts,
                cov = covFloat
        }
        call find_blocks_t.findBlocks {
            input:
                coverageGz = subsetCoverage.outputCoverageGz,
                table = fitModel.probabilityTable
        }
    }
    output {
        File bedsTarGz = select_first([makeNullOutput.bedsTarGz, findBlocks.bedsTarGz])
    }
}

task makeNullOutput{
    input{
        String prefix
        String suffix
        # runtime configurations
        Int memSize=2
        Int threadCount=2
        Int diskSize=8
        String dockerImage="quay.io/masri2019/hpp_base:latest"
        Int preemptible=2
    }
    command <<<
        mkdir ~{suffix}
        cd ~{suffix}
        touch ~{prefix}.~{suffix}.error.bed
        touch ~{prefix}.~{suffix}.duplicated.bed
        touch ~{prefix}.~{suffix}.haploid.bed
        touch ~{prefix}.~{suffix}.collapsed.bed
        cd ../

        tar -cf ~{prefix}.beds.~{suffix}.tar ~{suffix}
        gzip ~{prefix}.beds.~{suffix}.tar
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bedsTarGz = glob("*.tar.gz")[0]
    } 
}
