version 1.0

import "../tasks/long_read_aligner_scattered_PhaseHom.wdl" as long_read_aligner_scattered_t
import "./PHARAOH.wdl" as PHARAOH_t
import "../tasks/DeepPolisher.wdl" as deepPolisher_t
import "../tasks/applyPolish.wdl" as applyPolish_t
import "../tasks/parse_fastas.wdl" as parse_fastas_t
import "../tasks/filter_short_reads.wdl" as filter_short_reads_t

workflow hprc_DeepPolisher {
    meta {
        author: "Mira Mastoras"
        email: "mmastora@ucsc.edu"
        description: "Full HPRC polishing pipeline: aligning Hifi reads to raw diploid assembly, phasing homozygous regions with ONT UL, running DeepPolisher and applying polish to raw assemblies"
    }

    input {
        File Hap1RawFasta
        File Hap2RawFasta

        File DeepPolisherModelFilesTarGZ

        Array[File] ONTReads
        Array[File] HifiReads

        String DeepPolisherDocker="google/deepconsensus:polisher_v0.0.8_12122023"
        String sampleName
        Boolean useMargin=false

        # for minimap2, use k=19 for "map-hifi" and k=15 for "map-ont"
        # for winnowmap, use k=15 for preset "map-pb" and "map-ont"
        # default is minimap2

        String hifiAlignerToUse="minimap2"
        String ONTAlignerToUse="minimap2"
        String pafAligner="minimap2"
        String alignerHiFiPreset="map-hifi"
        String alignerONTPreset="map-ont"
        String alignerHiFiKmerSize="19"
        String alignerONTKmerSize="15"
    }

    ## parse input fasta files to obtain necessary forHap2s
    call parse_fastas_t.parseFastas as parseFastaStep {
        input:
            hap1Fasta=Hap1RawFasta,
            hap2Fasta=Hap2RawFasta,
            sampleName=sampleName
    }

    ## extract ONT reads > 100kb
    call filter_short_reads_t.FilterShortReads as extractONTUL100kb {
        input:
          readFiles=ONTReads,
          minReadLength=100000
    }

    ## Align all hifi reads to diploid assembly

    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignHifiToDiploid {
        input:
          assembly=parseFastaStep.dipRawFastaGz,
          readFiles=HifiReads,
          aligner=hifiAlignerToUse,
          preset=alignerHiFiPreset,
          kmerSize=alignerHiFiKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y -I8g",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="hifi.to.diploid.asm"
    }

    ## Align all ONT UL reads to Hap1 haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToHap1 {
        input:
          assembly=parseFastaStep.hap1RawFasta,
          readFiles=extractONTUL100kb.longReadFastqGzArray,
          aligner=ONTAlignerToUse,
          preset=alignerONTPreset,
          kmerSize=alignerONTKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="ONT.to.Hap1.asm"
    }

    ## Align all ONT UL reads to Hap2 haplotype
    call long_read_aligner_scattered_t.longReadAlignmentScattered as alignONTToHap2 {
        input:
          assembly=parseFastaStep.hap2RawFasta,
          readFiles=extractONTUL100kb.longReadFastqGzArray,
          aligner=ONTAlignerToUse,
          preset=alignerONTPreset,
          kmerSize=alignerONTKmerSize,
          sampleName=sampleName,
          options="--cs --eqx -L -Y",
          dockerImage="mobinasri/long_read_aligner:v0.3.3",
          sampleSuffix="ONT.to.Hap2.asm"
    }

    ## Phase reads in homozygous regions with UL, secphase marker mode in non-homoyzgous regions
    call PHARAOH_t.PHARAOH as PHARAOH {
        input:
          Hap1Fasta=parseFastaStep.hap1RawFasta,
          Hap1FastaIndex=parseFastaStep.hap1RawFastaIndex,
          Hap2Fasta=parseFastaStep.hap2RawFasta,
          Hap2FastaIndex=parseFastaStep.hap2RawFastaIndex,
          diploidFaGz=parseFastaStep.dipRawFastaGz,
          allHifiToDiploidBam=alignHifiToDiploid.bamFile,
          allHifiToDiploidBai=alignHifiToDiploid.baiFile,
          allONTToHap2Bam=alignONTToHap2.bamFile,
          allONTToHap1Bam=alignONTToHap1.bamFile,
          allONTToHap2Bai=alignONTToHap2.baiFile,
          allONTToHap1Bai=alignONTToHap1.baiFile,
          sampleName=sampleName,
          useMargin=useMargin,
          PharaohAligner=hifiAlignerToUse,
          PharaohHiFiPreset=alignerHiFiPreset,
          PharaohKmerSize=alignerHiFiKmerSize,
          pafAligner=pafAligner

    }

    ## Pass final phased hifi alignments to deepPolisher to produce polishing variants
    call deepPolisher_t.runDeepPolisher as DeepPolisher {
        input:
          Bam=PHARAOH.finalPhasedDipBam,
          Bai=PHARAOH.finalPhasedDipBai,
          Fasta=parseFastaStep.dipRawFastaGz,
          ModelFilesTarGZ=DeepPolisherModelFilesTarGZ,
          dockerImage=DeepPolisherDocker,
          sampleName=sampleName
    }

    ## Apply polishing variants to each haplotype
    call applyPolish_t.applyPolish as applyDPPolish_Hap1 {
        input:
          polishingVcf=DeepPolisher.PolisherVcf,
          asmRaw=parseFastaStep.hap1RawFasta,
          outPrefix=sampleName,
          HaplotypeLabel="Hap1"
    }

    call applyPolish_t.applyPolish as applyDPPolish_Hap2 {
        input:
          polishingVcf=DeepPolisher.PolisherVcf,
          asmRaw=parseFastaStep.hap2RawFasta,
          outPrefix=sampleName,
          HaplotypeLabel="Hap2"
    }

    output {
        File polishedAsmHap1=applyDPPolish_Hap1.asmPolished
        File polishedAsmHap2=applyDPPolish_Hap2.asmPolished
        File finalPhasedDipBam=PHARAOH.finalPhasedDipBam
        File finalPhasedDipBai=PHARAOH.finalPhasedDipBai
        File DeepPolisherVcf=DeepPolisher.PolisherVcf
        File DeepPolisherVcfUnfiltered=DeepPolisher.PolisherVcfNoFilter
        File variantBlocksBed=PHARAOH.variantBlocksBed
        File modifiedReadBlocksVariantsBed=PHARAOH.modifiedReadBlocksVariantsBed
        File DeepPolisherNoFiltersVcf=DeepPolisher.PolisherVcfNoFilter
    }
}
