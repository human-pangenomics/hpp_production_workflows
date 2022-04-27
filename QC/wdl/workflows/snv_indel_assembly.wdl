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

    call mergeVCF {
        input:
            VCF1     = filt_1.vcfOut,
            VCF2     = filt_2.vcfOut,
            happyVCF = happy_t.vcfOut,
            sample   = sample
    }

    call createVCFStats {
        input:
            inputVCF = mergeVCF.vcfOut
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

        File mergedVCF      = mergeVCF.vcfOut
        File mergedVCFIdx   = createVCFStats.vcfTBIOut
        File mergedVCFStats = createVCFStats.statsFileOut
    }
}

task mergeVCF {
    input{
        File VCF1
        File VCF2
        File happyVCF
        String sample

        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 128
        String dockerImage = "kishwars/t2t_polishing@sha256:418486a1e88c48555ad4f7158c0a9923762182e7c9cd883342ffe0a161d89de6" # v0.1


    }

    String outputVCF = "~{sample}_MERGED_SMALL_VARIANTS.vcf.gz"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace


        python3 vcf_merge_t2t.py \
            -v1 ~{VCF1} \
            -v2 ~{VCF2} \
            -hv ~{happyVCF} \
            -o ~{outputVCF}

    >>>
    output{
        File vcfOut = outputVCF
    }

    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task createVCFStats {
    input{
        File inputVCF

        Int memSizeGB = 8
        Int threadCount = 2
        Int diskSizeGB = 128
        String dockerImage = "kishwars/t2t_polishing@sha256:418486a1e88c48555ad4f7158c0a9923762182e7c9cd883342ffe0a161d89de6" # v0.1


    }

    String outputPrefix = basename(inputVCF, ".vcf.gz")
    String vcfTBI       = "~{outputPrefix}.vcf.gz.tbi"
    String statsFile    = "~{outputPrefix}.stats.txt"

    command <<<
        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        set -o xtrace

        ## Index VCF (not related to stats, just nice to have)
        tabix -p vcf ~{inputVCF}

        ## call bcftools stats
        bcftools stats ~{inputVCF} > ~{statsFile}
    >>>
    output{
        File vcfTBIOut    = vcfTBI
        File statsFileOut = statsFile
    }

    runtime{
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}