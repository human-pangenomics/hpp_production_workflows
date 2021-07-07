version 1.0

import "fit_model.wdl" as fitModel_t

workflow runFitModelMale {
    input {
        #mat
        File mat_autosome_nonCntr_Counts
        File mat_sex_nonCntr_Counts
        File mat_autosome_cntr_Counts
        File mat_sex_cntr_Counts
        #pat
        File pat_autosome_nonCntr_Counts
        File pat_sex_nonCntr_Counts
        File pat_autosome_cntr_Counts
        File pat_sex_cntr_Counts
    }
    scatter (matCounts in [mat_autosome_nonCntr_Counts, mat_sex_nonCntr_Counts, mat_autosome_cntr_Counts, mat_sex_cntr_Counts]){
        call fitModel_t.fitModel as matModel {
            input:
                counts = matCounts
        }
    }
    scatter (patCounts in [pat_autosome_nonCntr_Counts, pat_sex_nonCntr_Counts, pat_autosome_cntr_Counts, pat_sex_cntr_Counts]){
        call fitModel_t.fitModel as patModel {
            input:
                counts = patCounts
        }
    }
    output {
        #mat
        File mat_autosome_nonCntr_Table = matModel.probabilityTable[0]
        File mat_sex_nonCntr_Table = matModel.probabilityTable[1]
        File mat_autosome_cntr_Table = matModel.probabilityTable[2]
        File mat_sex_cntr_Table = matModel.probabilityTable[3]
        #pat
        File pat_autosome_nonCntr_Table = patModel.probabilityTable[0]
        File pat_sex_nonCntr_Table = patModel.probabilityTable[1]
        File pat_autosome_cntr_Table = patModel.probabilityTable[2]
        File pat_sex_cntr_Table = patModel.probabilityTable[3]
    }
}
