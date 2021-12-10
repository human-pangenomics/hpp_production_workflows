version 1.0

import "deep_variant.wdl" as dpv_t

workflow runDeepVariantNoSplit{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        String deepVariantModelType
        Int minMAPQ = 0
        String includeSecondary="False"
        String includeSupplementary="True"
    }
    call dpv_t.deepVariant{
        input:
            modelType = deepVariantModelType,
            assemblyFastaGz = assemblyFastaGz,
            bam = bam,
            includeSecondary = includeSecondary,
            includeSupplementary = includeSupplementary,
            minMAPQ = minMAPQ,
            threadCount=64,
            memSize=256,
            diskSize=512 
    }
    output{
        File vcfGz = deepVariant.vcfGz
    }
}

