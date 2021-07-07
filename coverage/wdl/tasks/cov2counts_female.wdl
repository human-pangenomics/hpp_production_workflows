version 1.0

import "cov2counts.wdl" as cov2counts_t

workflow runCov2CountsFemale {
    input {
        #mat
        File mat_autosome_nonCntr_CoverageGz
        File mat_autosome_cntr_CoverageGz 
        #pat
        File pat_autosome_nonCntr_CoverageGz
        File pat_autosome_cntr_CoverageGz
    }
    scatter (matCovGz in [mat_autosome_nonCntr_CoverageGz, mat_autosome_cntr_CoverageGz]){
        call cov2counts_t.cov2counts as matCounts {
            input:
                coverageGz = matCovGz
        }
    }
    scatter (patCovGz in [pat_autosome_nonCntr_CoverageGz, pat_autosome_cntr_CoverageGz]){
        call cov2counts_t.cov2counts as patCounts {
            input:
                coverageGz = patCovGz
        }
    }
    output {
        #mat
        File mat_autosome_nonCntr_Counts = matCounts.counts[0]
        File mat_autosome_cntr_Counts = matCounts.counts[1]
        #pat
        File pat_autosome_nonCntr_Counts = patCounts.counts[0]
        File pat_autosome_cntr_Counts = patCounts.counts[1]
    }
}
