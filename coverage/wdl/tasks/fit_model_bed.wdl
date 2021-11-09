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
            table = fitModel.probabilityTable,
            suffix = suffix 
    }
    output {
        File bedsTarGz = findBlocks.bedsTarGz
    }
}
