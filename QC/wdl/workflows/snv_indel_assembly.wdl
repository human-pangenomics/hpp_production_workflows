version 1.0

import "../tasks/deepvariant.wdl" as runDeepVariant
import "../tasks/pepperMarginDeepVariant.wdl" as runPepperMarginDeepVariant
import "../tasks/hapDotPy.wdl" as runHappy


workflow snv_indel_assembly {

    input {
        File assembly
        File assemblyIndex
        String sample
    }

    call runDeepVariant.DeepVariant as deepVariant_t {
        input:
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample
    }

    call runPepperMarginDeepVariant.pepperMarginDeepVariant as pmdv_t {
        input:
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample
    }    

    ## Filter DeepVariant calls
    call runPepperMarginDeepVariant.bcftoolsFilter as filt_1 {
        input:
            inputVCF      = deepVariant_t.vcfOut,
            excludeExpr   = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'",
            applyFilters  = "PASS"
    }

    ## Filter Pepper Margin DeepVariant calls
    call runPepperMarginDeepVariant.bcftoolsFilter as filt_2 {
        input:
            inputVCF      = pmdv_t.vcfOut,
            excludeExpr   = "'FORMAT/VAF<=0.5 | FORMAT/GQ<=30'",
            applyFilters  = "PASS",
            exludeTypes   = "indels"
    }

    call runHappy.hapDotPy as happy_t {
        input:    
            truthVCF      = filt_1.vcfOut,
            queryVCF      = filt_2.vcfOut,
            assembly      = assembly,
            assemblyIndex = assemblyIndex,
            sample        = sample        
    }

    output{
        File deepVariantVCF    = deepVariant_t.vcfOut
        File deepVariantVCFTBI = deepVariant_t.vcfIdxOut
        File deepVariantHTML   = deepVariant_t.visualReport

        File pepperMarginDeepVariantVCF     = pmdv_t.vcfOut
        File pepperMarginDeepVariantVCFIdx  = pmdv_t.vcfIdxOut
        File pepperMarginDeepVariantHTML    = pmdv_t.visualReport

        File DeepVariantVCF = filt_1.vcfOut
        File filtPmdvVCF    = filt_2.vcfOut

        File hapDotPyVCF    = happy_t.vcfOut
        File hapDotPyVCFIdx = happy_t.vcfIdxOut
        File hapDotPyTar    = happy_t.happyTar
    }
}