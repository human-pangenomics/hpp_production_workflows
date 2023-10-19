version 1.0

import "../../../assembly/wdl/tasks/len_filter_fasta.wdl" as filter_fasta_t
import "../../../assembly/wdl/tasks/break_into_contigs.wdl" as breakIntoContigs_t
import "../../../QC/wdl/workflows/short_qc.wdl" as shortQC_t
import "../../../QC/wdl/tasks/produce_fai.wdl" as produceFai_t
import "../../../QC/wdl/tasks/find_assembly_breakpoints.wdl" as findAssemblyBreakpoints_wf

workflow comparisonQC {

    input {
        File hap1Fasta # paternal
        File hap2Fasta # maternal
        String sampleName

        Int min_len  = 100000
        Boolean breakFasta = false

        ## QC Inputs
        File genesFasta
        File hs38Paf
        File matYak
        File patYak 
        File childYak

        ## findAssemblyBreakpoints Inputs
        File reference
        File annotationBed
        File annotationSD
        File annotationCENSAT

        Int memSizeQC     = 200
        Int diskSizeQC    = 256
        Int threadCountQC = 32
        Int preemptibleQC = 2
        String dockerQC   = "quay.io/masri2019/qc-stats:latest"
        String zonesQC    = "us-west2-a"
        
    }

    ## Find breakpoints for each assembly
    call findAssemblyBreakpoints_wf.findAssemblyBreakpoints as findAssemblyBreakpoints_hap1 {
        input:
            assembly=hap1Fasta,
            reference=reference,
            annotationBed=annotationBed,
            annotationSD=annotationSD,
            annotationCENSAT=annotationCENSAT
    }

    call findAssemblyBreakpoints_wf.findAssemblyBreakpoints as findAssemblyBreakpoints_hap2 {
        input:
            assembly=hap2Fasta,
            reference=reference,
            annotationBed=annotationBed,
            annotationSD=annotationSD,
            annotationCENSAT=annotationCENSAT
    }
    

    ## filter by size
    call filter_fasta_t.filter_fasta as filter_fasta_hap1 {
        input:
            inputFasta=hap1Fasta,
            sampleName=sampleName,
            outputFileTag="hap1_filt",
            min_size=min_len
    }


    call filter_fasta_t.filter_fasta as filter_fasta_hap2 {
        input:
            inputFasta=hap2Fasta,
            sampleName=sampleName,
            outputFileTag="hap2_filt",
            min_size=min_len
    }

    call shortQC_t.shortQC as shortQC_filt {
        input:
            patFasta    = filter_fasta_hap1.filteredFasta,
            matFasta    = filter_fasta_hap2.filteredFasta,

            patYak      = patYak,
            matYak      = matYak,
            childYak    = childYak,
            genesFasta  = genesFasta,
            hs38Paf     = hs38Paf,
            childID     = sampleName,
            memSize     = memSizeQC,
            threadCount = threadCountQC,
            dockerImage = dockerQC,
            preemptible = preemptibleQC,
            diskSize    = diskSizeQC,
            zones       = zonesQC
    }



    ## Create Fai for filtered fasta
    call produceFai_t.produceFai as produceFai_hap1 {
        input:
            assemblyGz=filter_fasta_hap1.filteredFasta
    }

    call produceFai_t.produceFai as produceFai_hap2 {
        input:
            assemblyGz=filter_fasta_hap2.filteredFasta
    }


    if (breakFasta) {
        ### break haplotype 1 on Ns ###
        call breakIntoContigs_t.breakIntoContigs as breakIntoContigs_hap1 {
            input:
                assemblyFaGz=filter_fasta_hap1.filteredFasta
        }

        ### break haplotype 2 on Ns ###
        call breakIntoContigs_t.breakIntoContigs as breakIntoContigs_hap2 {
            input:
                assemblyFaGz=filter_fasta_hap2.filteredFasta
        }

        call shortQC_t.shortQC as shortQC_filt_brk {
            input:
                patFasta    = breakIntoContigs_hap1.assemblyContigsFaGz,
                matFasta    = breakIntoContigs_hap2.assemblyContigsFaGz,

                patYak      = patYak,
                matYak      = matYak,
                childYak    = childYak,
                genesFasta  = genesFasta,
                hs38Paf     = hs38Paf,
                childID     = sampleName,
                memSize     = memSizeQC,
                threadCount = threadCountQC,
                dockerImage = dockerQC,
                preemptible = preemptibleQC,
                diskSize    = diskSizeQC,
                zones       = zonesQC
        }        

        ## Create Fai for filtered fasta
        call produceFai_t.produceFai as produceFai_hap1_brk {
            input:
                assemblyGz=breakIntoContigs_hap1.assemblyContigsFaGz
        }

        call produceFai_t.produceFai as produceFai_hap2_brk {
            input:
                assemblyGz=breakIntoContigs_hap2.assemblyContigsFaGz
        }

    }

	output {

        ## Outputs for breakpoint analysis
        File hap1_T2Tcontigs             = findAssemblyBreakpoints_hap1.T2Tcontigs
        File hap1_T2Tscaffolds           = findAssemblyBreakpoints_hap1.T2Tscaffolds
        File hap1_assembly_CHM13         = findAssemblyBreakpoints_hap1.assembly_CHM13
        File hap1_assemblyStatistics     = findAssemblyBreakpoints_hap1.assemblyStatistics
        File hap1_unifiedAssembly        = findAssemblyBreakpoints_hap1.unifiedAssembly
        File hap1_breakAnnotation_region = findAssemblyBreakpoints_hap1.breakAnnotation_region
        File hap1_breakAnnotation_SD     = findAssemblyBreakpoints_hap1.breakAnnotation_SD
        File hap1_breakAnnotation_CENSAT = findAssemblyBreakpoints_hap1.breakAnnotation_CENSAT
        File hap1_assembly_CHM13         = findAssemblyBreakpoints_hap1.assembly_CHM13 
        File hap1_filteredFlanksBed      = findAssemblyBreakpoints_hap1.filteredFlanksBed
        File hap1_bed_region             = findAssemblyBreakpoints_hap1.bed_region
        File hap1_bed_SD                 = findAssemblyBreakpoints_hap1.bed_SD
        File hap1_bed_CENSAT             = findAssemblyBreakpoints_hap1.bed_CENSAT
        File hap1_assemblyStatistics     = findAssemblyBreakpoints_hap1.assemblyStatistics

        File hap2_T2Tcontigs             = findAssemblyBreakpoints_hap2.T2Tcontigs
        File hap2_T2Tscaffolds           = findAssemblyBreakpoints_hap2.T2Tscaffolds
        File hap2_assembly_CHM13         = findAssemblyBreakpoints_hap2.assembly_CHM13
        File hap2_assemblyStatistics     = findAssemblyBreakpoints_hap2.assemblyStatistics
        File hap2_unifiedAssembly        = findAssemblyBreakpoints_hap2.unifiedAssembly
        File hap2_breakAnnotation_region = findAssemblyBreakpoints_hap2.breakAnnotation_region
        File hap2_breakAnnotation_SD     = findAssemblyBreakpoints_hap2.breakAnnotation_SD
        File hap2_breakAnnotation_CENSAT = findAssemblyBreakpoints_hap2.breakAnnotation_CENSAT
        File hap2_assembly_CHM13         = findAssemblyBreakpoints_hap2.assembly_CHM13 
        File hap2_filteredFlanksBed      = findAssemblyBreakpoints_hap2.filteredFlanksBed
        File hap2_bed_region             = findAssemblyBreakpoints_hap2.bed_region
        File hap2_bed_SD                 = findAssemblyBreakpoints_hap2.bed_SD
        File hap2_bed_CENSAT             = findAssemblyBreakpoints_hap2.bed_CENSAT
        File hap2_assemblyStatistics     = findAssemblyBreakpoints_hap2.assemblyStatistics


        ## Outputs for filtered assemblies
        File hap1FiltFasta     = filter_fasta_hap1.filteredFasta
        File hap2FiltFasta     = filter_fasta_hap2.filteredFasta

        File filtQCStats       = shortQC_filt.qcStatsText
        
        File hap1FiltFai       = produceFai_hap1.fai
        File hap2FiltFai       = produceFai_hap2.fai


        ## Outputs for broken assemblies
        File? hap1FiltBrkFasta = breakIntoContigs_hap1.assemblyContigsFaGz
        File? hap2FiltBrkFasta = breakIntoContigs_hap2.assemblyContigsFaGz

        File? filtBrkQCStats    = shortQC_filt_brk.qcStatsText

        File? hap1FiltBrkFai   = produceFai_hap1_brk.fai
        File? hap2FiltBrkFai   = produceFai_hap2_brk.fai
        

	}
}
