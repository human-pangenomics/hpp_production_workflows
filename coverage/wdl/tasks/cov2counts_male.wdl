version 1.0

import "cov2counts.wdl" as cov2counts_t

workflow runCov2CountsMale {
    input {
        #mat
        File mat_autosome_nonCntr_CoverageGz
        File mat_sex_nonCntr_CoverageGz
        File mat_autosome_cntr_CoverageGz 
        File mat_sex_cntr_CoverageGz
        #pat
        File pat_autosome_nonCntr_CoverageGz
        File pat_sex_nonCntr_CoverageGz
        File pat_autosome_cntr_CoverageGz
        File pat_sex_cntr_CoverageGz
    }
    scatter (matCovGz in [mat_autosome_nonCntr_CoverageGz, mat_sex_nonCntr_CoverageGz, mat_autosome_cntr_CoverageGz, mat_sex_cntr_CoverageGz]){
        call cov2counts_t.cov2counts as matCounts {
            input:
                coverageGz = matCovGz
        }
    }
    scatter (patCovGz in [pat_autosome_nonCntr_CoverageGz, pat_sex_nonCntr_CoverageGz, pat_autosome_cntr_CoverageGz, pat_sex_cntr_CoverageGz]){
        call cov2counts_t.cov2counts as patCounts {
            input:
                coverageGz = patCovGz
        }
    }
    output {
        #mat
        File mat_autosome_nonCntr_Counts = matCounts.counts[0]
        File mat_sex_nonCntr_Counts = matCounts.counts[1]
        File mat_autosome_cntr_Counts = matCounts.counts[2]
        File mat_sex_cntr_Counts = matCounts.counts[3]
        #pat
        File pat_autosome_nonCntr_Counts = patCounts.counts[0]
        File pat_sex_nonCntr_Counts = patCounts.counts[1]
        File pat_autosome_cntr_Counts = patCounts.counts[2]
        File pat_sex_cntr_Counts = patCounts.counts[3]
    }
}
