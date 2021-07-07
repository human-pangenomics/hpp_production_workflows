version 1.0

import "find_blocks.wdl" as findBlocks_t

workflow FindBlocksMale {
    input {
        #mat
        File mat_autosome_nonCntr_Table
        File mat_sex_nonCntr_Table
        File mat_autosome_cntr_Table
        File mat_sex_cntr_Table
        File mat_autosome_nonCntr_coverageGz
        File mat_sex_nonCntr_coverageGz
        File mat_autosome_cntr_coverageGz
        File mat_sex_cntr_coverageGz
        #pat
        File pat_autosome_nonCntr_Table
        File pat_sex_nonCntr_Table
        File pat_autosome_cntr_Table
        File pat_sex_cntr_Table
        File pat_autosome_nonCntr_coverageGz
        File pat_sex_nonCntr_coverageGz
        File pat_autosome_cntr_coverageGz
        File pat_sex_cntr_coverageGz
    }
    Array[File] matTables = [mat_autosome_nonCntr_Table, mat_sex_nonCntr_Table, mat_autosome_cntr_Table, mat_sex_cntr_Table]
    Array[File] matCoverages = [mat_autosome_nonCntr_coverageGz, mat_sex_nonCntr_coverageGz, mat_autosome_cntr_coverageGz, mat_sex_cntr_coverageGz]
    Array[File] patTables = [pat_autosome_nonCntr_Table, pat_sex_nonCntr_Table, pat_autosome_cntr_Table, pat_sex_cntr_Table]
    Array[File] patCoverages = [pat_autosome_nonCntr_coverageGz, pat_sex_nonCntr_coverageGz, pat_autosome_cntr_coverageGz, pat_sex_cntr_coverageGz] 
    scatter (matTableCov in zip(matTables, matCoverages)){
        call findBlocks_t.findBlocks as matBlocks {
            input:
                coverageGz = matTableCov.right,
                table =  matTableCov.left
        }
    }
    scatter (patTableCov in zip(patTables, patCoverages)){
        call findBlocks_t.findBlocks as patBlocks {
            input:
                coverageGz = patTableCov.right,
                table = patTableCov.left
        }
    }
    output {
        #mat
        Array[File] mat_autosome_nonCntr_Bed = matBlocks.bedFiles[0]
        Array[File] mat_sex_nonCntr_Bed = matBlocks.bedFiles[1]
        Array[File] mat_autosome_cntr_Bed = matBlocks.bedFiles[2]
        Array[File] mat_sex_cntr_Bed = matBlocks.bedFiles[3]
        #pat
        Array[File] pat_autosome_nonCntr_Bed = patBlocks.bedFiles[0]
        Array[File] pat_sex_nonCntr_Bed = patBlocks.bedFiles[1]
        Array[File] pat_autosome_cntr_Bed = patBlocks.bedFiles[2]
        Array[File] pat_sex_cntr_Bed = patBlocks.bedFiles[3]
    }
}
