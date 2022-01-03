version 1.0

import "../tasks/correct_bam.wdl" as correct_bam_t
import "../tasks/deep_variant.wdl" as deep_variant_t
import "../tasks/filter_alt_reads.wdl" as filter_alt_reads_t
import "../tasks/bam_coverage.wdl" as bam_coverage_t
import "../tasks/pepper_margin_deep_variant.wdl" as pmdv_t

workflow runCoverageAnalysisStep1{
    input {
        File bam
        File assemblyFastaGz
        File? phasingLogText
        Int minReadLength 
        Int minAlignmentLength
        String deepVariantModelType
        String pepperModelType
        String variantCaller = "dv"
    }
    
    ## Correct the bam file by swapping pri/sec tags for the wrongly phased reads
    call correct_bam_t.correctBam {
        input:        
            bam = bam,
            phasingLogText = phasingLogText,
            suffix = "corrected",
            options = "--primaryOnly --minReadLen ${minReadLength} --minAlignment ${minAlignmentLength}",
            diskSize = ceil(size(bam, "GB")) * 2 + 64
    }
    
    if ("${variantCaller}" == "dv") { 
        ## Call variants to be used for finding the reads with alternative alleles
        call deep_variant_t.runVariantCalling as dpv {
            input:
                deepVariantModelType = deepVariantModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = correctBam.correctedBam,
                bamIndex = correctBam.correctedBamIndex,
                minMAPQ = 0,
                includeSecondary="False",
                includeSupplementary="True"
        }
    }
    if ("${variantCaller}" == "pmdv") {
        ## Call variants to be used for finding the reads with alternative alleles
        call pmdv_t.runPepperMarginDeepVariant as pmdv {
            input:
                pmdvModelType = deepVariantModelType,
                assemblyFastaGz = assemblyFastaGz,
                bam = correctBam.correctedBam,
                bamIndex = correctBam.correctedBamIndex,
                minMAPQ = 0,
                includeSupplementary="True"
        }
    }
   
    File vcfGz1 = select_first([dpv.vcfGz, pmdv.vcfGz])     

    ## Filter the reads with alternative alleles
    call filter_alt_reads_t.filterAltReads {
        input:
            vcf = vcfGz1,
            bam = correctBam.correctedBam,
            diskSize = ceil(size(correctBam.correctedBam, "GB")) * 2 + 64
    }
    
    ## Calculate coverage for the corrected bam file (without filtering)
    call bam_coverage_t.bamCoverage as bam2cov_corrected{
        input:
            bam = correctBam.correctedBam,
            minMAPQ = 0,
            assemblyFastaGz = assemblyFastaGz,
            diskSize = ceil(size(correctBam.correctedBam, "GB")) + 256
    }
    
    ## Calculate coverage for the corrected bam file in which the 
    ## reads with alternative alleles are removed
    call bam_coverage_t.bamCoverage as bam2cov_altRemoved{
        input:
            bam = filterAltReads.filteredBam,
            minMAPQ = 0,
            assemblyFastaGz = assemblyFastaGz,
            diskSize = ceil(size(filterAltReads.filteredBam, "GB")) + 256
    }

    ## Calculate coverage of reads with high mapqs (20<) for the 
    ## corrected bam file (without filtering)
    ##
    ## This coverage will be used for checking the false duplications
    ## in the 2nd step of the coverage analysis
    call bam_coverage_t.bamCoverage as bam2cov_corrected_highMapq{
        input:
            bam = correctBam.correctedBam,
            minMAPQ = 20,
            assemblyFastaGz = assemblyFastaGz,
            diskSize = ceil(size(correctBam.correctedBam, "GB")) + 256
    }

    ## Calculate coverage of reads with high mapqs (20<) for the 
    ## corrected bam file in which the
    ## reads with alternative alleles are removed
    ##
    ## This coverage will be used for checking the false duplications
    ## in the 2nd step of the coverage analysis
    call bam_coverage_t.bamCoverage as bam2cov_altRemoved_highMapq{
        input:
            bam = filterAltReads.filteredBam,
            minMAPQ = 20,
            assemblyFastaGz = assemblyFastaGz,
            diskSize = ceil(size(filterAltReads.filteredBam, "GB")) + 256
    }
    output {
        File vcfGz = vcfGz1
        File altBam = filterAltReads.altBam
        File altBai = filterAltReads.altBamIndex
        Float meanCorrectedCoverageFloat = bam2cov_corrected.coverageMean
        Float meanAltRemovedCoverageFloat = bam2cov_altRemoved.coverageMean
        File correctedCovGz = bam2cov_corrected.coverageGz
        File altRemovedCovGz = bam2cov_altRemoved.coverageGz
        File correctedHighMapqCovGz = bam2cov_corrected_highMapq.coverageGz
        File altRemovedHighMapqCovGz = bam2cov_altRemoved_highMapq.coverageGz
    }
}
