version 1.0

import "bedtools.wdl" as bedtools_t

workflow produceSexBed {
    input {
        File chrX_nonCntr_nonPAR_bed
        File chrY_nonPAR_bed
        String assemblyFastaGz
        String suffix
    }
    call bedtools_t.union {
       input:
          bedFiles = [chrX_nonCntr_nonPAR_bed, chrY_nonPAR_bed],
          outputPrefix = basename("${assemblyFastaGz}", ".fa.gz") + ".${suffix}"
    }
    output {
       File sex_nonCntr_nonPAR_bed = union.unionBed
    }
}
