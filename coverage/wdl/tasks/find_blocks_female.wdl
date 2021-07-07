version 1.0

import "find_blocks.wdl" as findBlocks_t

workflow runFitModelFemale {
    input {
        #mat
        File mat_autosome_nonCntr_Table
        File mat_autosome_cntr_Table
        File mat_autosome_nonCntr_CoverageGz
        File mat_autosome_cntr_CoverageGz
        #pat
        File pat_autosome_nonCntr_Table
        File pat_autosome_cntr_Table
        File pat_autosome_nonCntr_CoverageGz
        File pat_autosome_cntr_CoverageGz
    }
    Array[File] matTables = [mat_autosome_nonCntr_Table, mat_autosome_cntr_Table]
    Array[File] patTables = [pat_autosome_nonCntr_Table, pat_autosome_cntr_Table]
    Array[File] matCoverages = [mat_autosome_nonCntr_CoverageGz, mat_autosome_cntr_CoverageGz]
    Array[File] patCoverages = [pat_autosome_nonCntr_CoverageGz, pat_autosome_cntr_CoverageGz]
    
   scatter (matTableCov in zip(matTables, matCoverages)){
        call findBlocks_t.findBlocks as matBlocks {
            input:
                table = matTableCov.left,
                coverageGz = matTableCov.right
        }
    }
    scatter (patTableCov in zip(patTables, patCoverages)){
        call findBlocks_t.findBlocks as patBlocks {
            input:
                table = patTableCov.left,
                coverageGz = patTableCov.right
        }
    }
    output {
        #mat
        Array[File] mat_autosome_nonCntr_Bed = matBlocks.bedFiles[0]
        Array[File] mat_autosome_cntr_Bed = matBlocks.bedFiles[1]
        #pat
        Array[File] pat_autosome_nonCntr_Bed = patBlocks.bedFiles[0]
        Array[File] pat_autosome_cntr_Bed = patBlocks.bedFiles[1]
    }
}
