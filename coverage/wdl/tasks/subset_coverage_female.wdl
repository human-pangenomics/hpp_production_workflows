version 1.0

import "subset_coverage.wdl" as scov_t

workflow runSubsetCoverageMale {
    input {
        File autosome_nonCntr_bed
        File autosome_cntr_bed
        File coverageGz
    }
    call scov_t.subsetCoverage as sub_1 {
        input:
            coverageGz = coverageGz,
            blocksBed = autosome_nonCntr_bed,
            suffix = "autosome_nonCntr"
    }
    call scov_t.subsetCoverage as sub_2 {
        input:
            coverageGz = coverageGz,
            blocksBed = autosome_cntr_bed,
            suffix = "autosome_cntr"
    }
    output {
       File autosome_nonCntr_cov_gz = sub_1.outputCoverageGz
       File autosome_cntr_filled_cov_gz = sub_2.outputCoverageGz
    }
}
