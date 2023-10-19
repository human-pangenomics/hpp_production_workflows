version 1.0

import "../../../assembly/wdl/tasks/len_filter_fasta.wdl" as filter_fasta_t
import "../../../assembly/wdl/tasks/break_into_contigs.wdl" as breakIntoContigs_t
import "../../../QC/wdl/workflows/short_diploid_qc.wdl" as shortQC_t
import "../../../QC/wdl/tasks/produce_fai.wdl" as produceFai_t
import "../../../QC/wdl/tasks/find_assembly_breakpoints.wdl" as findAssemblyBreakpoints_wf

workflow comparisonQC {

    input {
        File asmFasta  # diploid
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
        String dockerQC   = "mobinasri/hpp_qc_stats:latest"
        String zonesQC    = "us-west2-a"
        
    }

    ## Find breakpoints
    call findAssemblyBreakpoints_wf.findAssemblyBreakpoints as findAssemblyBreakpoints_t {
        input:
            assembly=asmFasta,
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

    call filter_fasta_t.filter_fasta as filter_fasta_dip {
        input:
            inputFasta=asmFasta,
            sampleName=sampleName,
            outputFileTag="dip_filt",
            min_size=min_len
    }

    call shortQC_t.shortDiploidQC as shortQC_filt {
        input:
            patFastaGz  = filter_fasta_hap1.filteredFasta,
            matFastaGz  = filter_fasta_hap2.filteredFasta,
            dipFastaGz  = filter_fasta_dip.filteredFasta,

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

        ### break diploid on Ns ###
        call breakIntoContigs_t.breakIntoContigs as breakIntoContigs_dip {
            input:
                assemblyFaGz=filter_fasta_dip.filteredFasta
        }

        call shortQC_t.shortDiploidQC as shortQC_filt_brk {
            input:
                patFastaGz  = breakIntoContigs_hap1.assemblyContigsFaGz,
                matFastaGz  = breakIntoContigs_hap2.assemblyContigsFaGz,
                dipFastaGz  = breakIntoContigs_dip.assemblyContigsFaGz,

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

        ## Create Fai for filtered, broken fasta
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
        File dip_T2Tcontigs             = findAssemblyBreakpoints_t.T2Tcontigs
        File dip_T2Tscaffolds           = findAssemblyBreakpoints_t.T2Tscaffolds
        File dip_assembly_CHM13         = findAssemblyBreakpoints_t.assembly_CHM13
        File dip_assemblyStatistics     = findAssemblyBreakpoints_t.assemblyStatistics
        File dip_unifiedAssembly        = findAssemblyBreakpoints_t.unifiedAssembly
        File dip_breakAnnotation_region = findAssemblyBreakpoints_t.breakAnnotation_region
        File dip_breakAnnotation_SD     = findAssemblyBreakpoints_t.breakAnnotation_SD
        File dip_breakAnnotation_CENSAT = findAssemblyBreakpoints_t.breakAnnotation_CENSAT
        File dip_assembly_CHM13         = findAssemblyBreakpoints_t.assembly_CHM13 
        File dip_filteredFlanksBed      = findAssemblyBreakpoints_t.filteredFlanksBed
        File dip_bed_region             = findAssemblyBreakpoints_t.bed_region
        File dip_bed_SD                 = findAssemblyBreakpoints_t.bed_SD
        File dip_bed_CENSAT             = findAssemblyBreakpoints_t.bed_CENSAT
        File dip_assemblyStatistics     = findAssemblyBreakpoints_t.assemblyStatistics


        ## Outputs for filtered assemblies
        File hap1FiltFasta     = filter_fasta_hap1.filteredFasta
        File hap2FiltFasta     = filter_fasta_hap2.filteredFasta
        File dipFiltFasta      = filter_fasta_dip.filteredFasta

        File filtQCStats       = shortQC_filt.qcStatsText
        
        File hap1FiltFai       = produceFai_hap1.fai
        File hap2FiltFai       = produceFai_hap2.fai


        ## Outputs for broken assemblies
        File? hap1FiltBrkFasta = breakIntoContigs_hap1.assemblyContigsFaGz
        File? hap2FiltBrkFasta = breakIntoContigs_hap2.assemblyContigsFaGz
        File? dipFiltBrkFasta  = breakIntoContigs_dip.assemblyContigsFaGz
 
        File? filtBrkQCStats    = shortQC_filt_brk.qcStatsText

        File? hap1FiltBrkFai   = produceFai_hap1_brk.fai
        File? hap2FiltBrkFai   = produceFai_hap2_brk.fai
        

	}
}
