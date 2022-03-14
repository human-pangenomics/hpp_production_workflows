version 1.0

import "deep_variant.wdl" as dv_t
import "../../../asset/wdl/tasks/variant_calling.wdl" as var_t
import "pepper_margin_deep_variant.wdl" as pmdv_t

workflow runPepperMarginDeepVariantSplit{
    input {
        File assemblyFastaGz
        File bam
        File bamIndex
        ## Model types can be:
        ##     "ont_r9_guppy5_sup"
        ##     "ont_r10_q20"
        ##     "hifi"
        String pmdvModelType
        Int minMAPQ = 0
        String includeSupplementary="True"
        Int numberOfCallerNodes=16
        Int nodeThreadCount=8 
    }
    call pmdv_t.removeMultiplePrimary{
        input:
            bam = bam,
            diskSize = 2 * ceil(size(bam, "GB")) + 64
    }
    call dv_t.splitBamContigWise {
        input:
            assemblyFastaGz = assemblyFastaGz,
            bam = removeMultiplePrimary.correctedBam,
            bamIndex = removeMultiplePrimary.correctedBai,
            splitNumber = numberOfCallerNodes,
            threadCount = numberOfCallerNodes,
            diskSize = 2 * ceil(size(bam, "GB")) + 64
    }
    scatter (part in zip(splitBamContigWise.splitBams, splitBamContigWise.splitBeds)) {
        call pmdv_t.pmdv{
            input:
                modelType = pmdvModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = part.left,
                includeSupplementary = includeSupplementary,
                minMAPQ = minMAPQ,
                threadCount = nodeThreadCount,
                memSize = 16,
                diskSize= 2 * ceil(size(part.left, "GB")) + 64
        }
    }
    call var_t.mergeVcf{
        input:
            vcfGzFiles = pmdv.vcfGz,
            outputName = basename("${bam}", ".bam")
    }
    output{
        File vcfGz = mergeVcf.vcfGz
        Array[File] intermediateTars = pmdv.intermediateTar
    }
}
