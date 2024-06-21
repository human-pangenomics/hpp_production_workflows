version 1.0

import "../tasks/long_read_aligner.wdl" as long_read_aligner_t
import "../tasks/find_homozygous_regions.wdl" as findHomozygousRegions_t
import "../tasks/subBamByBed.wdl" as subBamByBed_t
import "../tasks/extract_reads.wdl" as extract_reads_t
import "../tasks/correct_bam.wdl" as correct_bam_t
import "../tasks/deepvariant.wdl" as deepvariant_t
import "../tasks/pepperMarginDeepVariant.wdl" as pmdv_t
import "../tasks/whatsHapPhase.wdl" as whatshap_phase_t
import "../tasks/marginPhase.wdl" as margin_phase_t
import "../tasks/long_read_aligner_scattered_PhaseHom.wdl" as long_read_aligner_scattered_t
import "../tasks/secphase.wdl" as secphase_t
import "../tasks/concatVcf.wdl" as concatVcf_t
import "../tasks/get_mapq_table.wdl" as get_mapq_t


workflow PHARAOH{
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "PHAsing Reads in Areas Of Homozygosity using ONT UL data"
    }
    input {
        File Hap1Fasta
        File Hap1FastaIndex
        File Hap2Fasta
        File Hap2FastaIndex

        File diploidFaGz

        File allHifiToDiploidBam
        File allHifiToDiploidBai

        File allONTToHap2Bam
        File allONTToHap1Bam
        File allONTToHap2Bai
        File allONTToHap1Bai

        String sampleName

        # PHARAOH defaults to using WhatsHap for phasing. To use Margin instead, set to true
        Boolean useMargin = false

        # option to pass in separate config to margin
        File? marginConfig

        # for minimap2, use k=19 for "map-hifi" and k=15 for "map-ont"
        # for winnowmap, use k=15 for preset "map-pb" and "map-ont"
        # default is minimap2

        String PharaohAligner="minimap2"
        String PharaohKmerSize=19
        String PharaohHiFiPreset="map-hifi"
        String pafAligner="minimap2"

        String minWindowSizeBp=20000
        String extendBp=50000
        String hifiAlignmentOptions="--cs --eqx -Y -L"
    }

    ## Align Hap2 to Hap1 assembly
    call long_read_aligner_t.alignmentPaf as alignmentPaf{
        input:
            aligner=pafAligner,
            preset="asm5",
            options="-L --eqx --cs -c",
            readFastq_or_queryAssembly=Hap2Fasta,
            refAssembly=Hap1Fasta,
            suffix="hap2_to_hap1",
            diskSize=512,
            threadCount=32,
            kmerSize=19,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    ## Get Homozygous regions between the two haplotypes
    call findHomozygousRegions_t.FindHomozygousRegions as findHomozygousRegions{
        input:
            pafFile=alignmentPaf.pafFile,
            minWindowSizeBp=minWindowSizeBp,
            extendBp=extendBp,
            outPrefix=sampleName
    }

    ## subset diploid bamfile to homozygous regions
    call subBamByBed_t.SubBamByBed as subDipBamByHomozygous{
        input:
            Bam=allHifiToDiploidBam,
            Bai=allHifiToDiploidBai,
            Bed=findHomozygousRegions.extendedBed
    }

    ## Align all homozygous reads to both haplotypes separately
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignAllToHap1Scattered{
        input:
            readFiles=[subDipBamByHomozygous.subBam],
            assembly=Hap1Fasta,
            aligner=PharaohAligner,
            preset=PharaohHiFiPreset,
            sampleName=sampleName,
            kmerSize=PharaohKmerSize,
            sampleSuffix="all2hap1",
            options=hifiAlignmentOptions,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignAllToHap2Scattered{
        input:
            readFiles=[subDipBamByHomozygous.subBam],
            assembly=Hap2Fasta,
            aligner=PharaohAligner,
            preset=PharaohHiFiPreset,
            kmerSize=PharaohKmerSize,
            sampleName=sampleName,
            sampleSuffix="all2hap2",
            options=hifiAlignmentOptions,
            dockerImage="mobinasri/long_read_aligner:v0.3.3"
    }

    ## remove reads with de > 0.02
    call correct_bam_t.correctBam as correctBamMaxDivergenceHap1 {
        input:
            Bam=alignAllToHap1Scattered.bamFile,
            options="--maxDiv 0.02",
            suffix="maxDiv.02",
            dockerImage="mobinasri/secphase:dev-v0.2.0-hom"

    }

    call correct_bam_t.correctBam as correctBamMaxDivergenceHap2 {
        input:
            Bam=alignAllToHap2Scattered.bamFile,
            options="--maxDiv 0.02",
            suffix="maxDiv.02",
            dockerImage="mobinasri/secphase:dev-v0.2.0-hom"

    }

    ## call variants on each bam
    call deepvariant_t.DeepVariant as DeepVariantHap1{
        input:
            inputReads=correctBamMaxDivergenceHap1.correctedBam,
            inputReadsIdx=correctBamMaxDivergenceHap1.correctedBamIndex,
            assembly=Hap1Fasta,
            assemblyIndex=Hap1FastaIndex,
            sample=sampleName,
            modelType = "PACBIO"
    }
    call deepvariant_t.DeepVariant as DeepVariantHap2{
        input:
            inputReads=correctBamMaxDivergenceHap2.correctedBam,
            inputReadsIdx=correctBamMaxDivergenceHap2.correctedBamIndex,
            assembly=Hap2Fasta,
            assemblyIndex=Hap2FastaIndex,
            sample=sampleName,
            modelType = "PACBIO"
    }

    ## filter variants by GQ
    call pmdv_t.bcftoolsFilter as FilterDVHap1{
        input:
          inputVCF=DeepVariantHap1.vcfOut,
          excludeExpr="'FORMAT/GQ<=10'",
          applyFilters=""
    }
    call pmdv_t.bcftoolsFilter as FilterDVHap2{
        input:
          inputVCF=DeepVariantHap2.vcfOut,
          excludeExpr="'FORMAT/GQ<=10'",
          applyFilters=""
    }

    ## Phase variants with UL reads using Margin or WhatsHap

    if (useMargin==false) {
        call whatshap_phase_t.WhatsHapPhase as WhatsHapPhaseHap2 {
            input:
              vcfFile=FilterDVHap2.vcfOut,
              vcfFileIdx=FilterDVHap2.vcfOutIdx,
              refFile=Hap2Fasta,
              refFileIdx=Hap2FastaIndex,
              bamFile=allONTToHap2Bam,
              bamFileIdx=allONTToHap2Bai,
              outPrefix="phased_Vcf_UL_Hap2"
        }

        call whatshap_phase_t.WhatsHapPhase as WhatsHapPhaseHap1 {
            input:
              vcfFile=FilterDVHap1.vcfOut,
              vcfFileIdx=FilterDVHap1.vcfOutIdx,
              refFile=Hap1Fasta,
              refFileIdx=Hap1FastaIndex,
              bamFile=allONTToHap1Bam,
              bamFileIdx=allONTToHap1Bai,
              outPrefix="phased_Vcf_UL_Hap1"
        }
    }

    if (useMargin==true) {
        call margin_phase_t.marginPhase as marginPhaseHap1 {
            input:
              vcfFile=FilterDVHap1.vcfOut,
              vcfFileIdx=FilterDVHap1.vcfOutIdx,
              refFile=Hap1Fasta,
              refFileIdx=Hap1FastaIndex,
              bamFile=allONTToHap1Bam,
              bamFileIdx=allONTToHap1Bai,
              outPrefix="phased_Vcf_UL_Hap1",
              HifiOrONT="ONT",
              configFile=marginConfig
        }

        call margin_phase_t.marginPhase as marginPhaseHap2 {
            input:
              vcfFile=FilterDVHap2.vcfOut,
              vcfFileIdx=FilterDVHap2.vcfOutIdx,
              refFile=Hap2Fasta,
              refFileIdx=Hap2FastaIndex,
              bamFile=allONTToHap2Bam,
              bamFileIdx=allONTToHap2Bai,
              outPrefix="phased_Vcf_UL_Hap2",
              HifiOrONT="ONT",
              configFile=marginConfig
        }
    }

    call concatVcf_t.bcftoolsConcat as bcftoolsConcat {
        input:
          vcf1 = select_first([marginPhaseHap1.phasedVcf, WhatsHapPhaseHap1.phasedVcf]) ,
          vcf2 = select_first([marginPhaseHap2.phasedVcf, WhatsHapPhaseHap2.phasedVcf])
    }
    call secphase_t.runSecPhase as runSecPhase {
        input:
          inputBam=allHifiToDiploidBam,
          diploidAssemblyFastaGz=diploidFaGz,
          phasedVcf=bcftoolsConcat.vcfOut,
          variantBed=findHomozygousRegions.bed
    }

    call get_mapq_t.runGetMapQTable as getMapQ {
        input:
          allHifiToHap2Bam=correctBamMaxDivergenceHap2.correctedBam,
          allHifiToHap2Bai=correctBamMaxDivergenceHap2.correctedBamIndex,
          allHifiToHap1Bam=correctBamMaxDivergenceHap1.correctedBam,
          allHifiToHap1Bai=correctBamMaxDivergenceHap1.correctedBamIndex,
          secPhaseBed=runSecPhase.variantBlocksBed
    }

    call correct_bam_t.correctBam as correctBamSecPhase {
        input:
          Bam=allHifiToDiploidBam,
          mapqTableText=getMapQ.mapqTable,
          phasingLogText=runSecPhase.outLog,
          suffix="PHARAOH",
          options="--maxDiv 0.002"
    }

    output {
        File asm2asmPaf=alignmentPaf.pafFile

        File phasedVcf=bcftoolsConcat.vcfOut

        File homExtendedbed=findHomozygousRegions.extendedBed
        File homBed=findHomozygousRegions.bed

        File dipBamHomozygous=subDipBamByHomozygous.subBam

        File allHomToHap1Bam=alignAllToHap1Scattered.bamFile
        File allHomToHap2Bam=alignAllToHap2Scattered.bamFile

        File allHomToHap1BamMaxDiv=correctBamMaxDivergenceHap1.correctedBam
        File allHomToHap2BamMaxDiv=correctBamMaxDivergenceHap2.correctedBam

        File deepVariantHap1=DeepVariantHap1.vcfOut
        File deepVariantHap2=DeepVariantHap2.vcfOut

        File deepVariantHap1Filt=FilterDVHap1.vcfOut
        File deepVariantHap2Filt=FilterDVHap2.vcfOut

        File secphaseOutLog=runSecPhase.outLog
        File variantBlocksBed=runSecPhase.variantBlocksBed
        File modifiedReadBlocksVariantsBed=runSecPhase.modifiedReadBlocksVariantsBed

        File finalPhasedDipBam=correctBamSecPhase.correctedBam
        File finalPhasedDipBai=correctBamSecPhase.correctedBamIndex


    }
}
