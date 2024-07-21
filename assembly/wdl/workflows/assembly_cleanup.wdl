version 1.0

import "../tasks/mitoHiFi.wdl" as mitoHiFi_wf
import "../../../QC/wdl/workflows/FCS.wdl" as fcs_wf
import "../../../QC/wdl/tasks/findMitoContigs.wdl" as findMito_wf
import "../../../QC/wdl/tasks/dropFastaContigs.wdl" as dropContigs_wf
import "https://raw.githubusercontent.com/biomonika/HPP/14474a7ac061e6f1aa7faa1d9ea19217a3efe8d1/assembly/wdl/workflows/evaluateHumanAssembly.wdl" as evaluateHumanAssembly_wf
import "../tasks/len_filter_fasta.wdl" as filter_fasta_t
import "../tasks/renameContigsAddMT.wdl" as renameContigsAddMito

workflow assembly_cleanup_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Workflow to take an assembly from the assembler to Genbank upload."
    }
    
    input {
        File hap1_fasta_gz
        File hap2_fasta_gz 
        String sample_id

        ## for mitoHiFi & findMitoContigs
        File related_mito_fasta

        ## for mitoHiFi
        File related_mito_genbank
        Array[File] hifi_reads
        File? referenceFasta
        Int fileExtractionDiskSizeGB = 512

        ## for findMitoContigs
        File ncbi_mito_blast_db

        ## for NCBI's foreign contamination screen
        File blast_div
        File GXI
        File GXS
        File manifest
        File metaJSON
        File seq_info
        File taxa

        ## For length filtering fastas
        Int min_sequence_len = 50000
    }


    ## mitoHiFi to assemble mitochondrial contig
    call mitoHiFi_wf.mitoHifiWorkflow as assemble_mito {
        input:
            sample_id                = sample_id,
            related_mito_fasta       = related_mito_fasta,
            related_mito_genbank     = related_mito_genbank,
            hifi_reads               = hifi_reads,
            filter_hifi              = true,
            referenceFasta           = referenceFasta,
            fileExtractionDiskSizeGB = fileExtractionDiskSizeGB
    }

    ## NCBI's Foreign Contamination Screen: Hap 1
    call fcs_wf.RunFCS as contam_screen_hap1 {
        input:
            assembly   = hap1_fasta_gz,
            asm_name   = "~{sample_id}_hap1",

            blast_div  = blast_div,
            GXI        = GXI,
            GXS        = GXS,
            manifest   = manifest,
            metaJSON   = metaJSON,
            seq_info   = seq_info,
            taxa       = taxa
    }

    ## NCBI's Foreign Contamination Screen: Hap 2
    call fcs_wf.RunFCS as contam_screen_hap2 {
        input:
            assembly   = hap2_fasta_gz,
            asm_name   = "~{sample_id}_hap2",

            blast_div  = blast_div,
            GXI        = GXI,
            GXS        = GXS,
            manifest   = manifest,
            metaJSON   = metaJSON,
            seq_info   = seq_info,
            taxa       = taxa
    } 


    call findMito_wf.findMitoContigs as findMitosInHap1 {
        input:
            sample_id          = "~{sample_id}",
            haplotype          = "hap1",
            inputFastaGZ       = contam_screen_hap1.output_fasta,
            mito_ref           = related_mito_fasta,
            ncbi_mito_blast_db = ncbi_mito_blast_db
    }

    call findMito_wf.findMitoContigs as findMitosInHap2 {
        input:
            sample_id          = "~{sample_id}",
            haplotype          = "hap2",
            inputFastaGZ       = contam_screen_hap2.output_fasta,
            mito_ref           = related_mito_fasta,
            ncbi_mito_blast_db = ncbi_mito_blast_db
    }

    ## remove mitos from hap1 (will add back in MitoHiFi version later)
    call dropContigs_wf.dropContigs as dropHap1Mito {
        input:
            sampleName        = "~{sample_id}",
            outputFileTag     = "hap1_mito_stripped", 
            inputFasta        = contam_screen_hap1.output_fasta,
            contigsToDrop     = findMitosInHap1.mitoContigs
    }  

    ## remove mitos from hap2 (will add back in MitoHiFi version later)
    call dropContigs_wf.dropContigs as dropHap2Mito {
        input:
            sampleName        = "~{sample_id}",
            outputFileTag     = "hap2_mito_stripped", 
            inputFasta        = contam_screen_hap2.output_fasta,
            contigsToDrop     = findMitosInHap2.mitoContigs
    }  


    call evaluateHumanAssembly_wf.evaluateHumanAssembly as evaluateHap1 {
        input:
            assembly          = dropHap1Mito.FinalAssembly
    }

    call evaluateHumanAssembly_wf.evaluateHumanAssembly as evaluateHap2 {
        input:
            assembly          = dropHap2Mito.FinalAssembly
    }

    ## filter by size
    call filter_fasta_t.filter_fasta as filter_fasta_hap1 {
        input:
            inputFasta=evaluateHap1.unifiedAssembly,
            sampleName=sample_id,
            outputFileTag="hap1_filt",
            min_size=min_sequence_len
    }   

    call filter_fasta_t.filter_fasta as filter_fasta_hap2 {
        input:
            inputFasta=evaluateHap2.unifiedAssembly,
            sampleName=sample_id,
            outputFileTag="hap2_filt",
            min_size=min_sequence_len
    }   


    ## rename contigs for Hap1
    call renameContigsAddMito.renameContigsAddMT as renameHeadersHap1 {
        input:
            sampleName     = sample_id,
            outputFileTag  = "hap1_for_genbank",
            haplotype      = 1,
            inputFastaGZ   = filter_fasta_hap1.filteredFasta,
            t2t_sequences  = evaluateHap1.T2Tscaffolds
    }

    ## rename contigs and add mito for Hap2
    call renameContigsAddMito.renameContigsAddMT as renameHeadersAddMitoHap2 {
        input:
            sampleName     = sample_id,
            outputFileTag  = "hap2_for_genbank",
            haplotype      = 2,
            inputFastaGZ   = filter_fasta_hap2.filteredFasta,
            t2t_sequences  = evaluateHap2.T2Tscaffolds,
            mitoAssembly   = assemble_mito.mitoHiFi_assembly
    }



    output {

        ## Final Assembly: Ready for Genbank
        File hap1_output_fasta_gz           = renameHeadersHap1.FinalAssembly
        File hap2_output_fasta_gz           = renameHeadersAddMitoHap2.FinalAssembly


        ## MitoHiFi outputs
        File mitoHiFi_assembly               = assemble_mito.mitoHiFi_assembly
        File mitoHiFi_stats                  = assemble_mito.mitoHiFi_stats
        File mitoHiFi_output_tar             = assemble_mito.mitoHiFi_output_tar
        
        ## paf of read alignments to concatenated mito (for troublshooting cutoffs)
        Array[File?] read_to_concatmito_paf  = assemble_mito.read_to_concatmito_paf

        ## mito assembly evaluations
        File mitoHiFi_eval_bam               = assemble_mito.mitoHiFi_eval_bam
        File mitoHiFi_eval_vcf               = assemble_mito.mitoHiFi_eval_vcf


        ## FCS outputs; hap1
        File fcs_hap1_int_clean_fa   = contam_screen_hap1.intermediate_clean_fa
        File fcs_hap1_contamFasta    = contam_screen_hap1.contamFasta
        File fcs_hap1_gx_report      = contam_screen_hap1.fcs_gx_report
        File fcs_hap1_tax_report     = contam_screen_hap1.fcs_taxonomy_report
        File fcs_hap1_adapter_report = contam_screen_hap1.adapter_Report
        File fcs_hap1_output_fasta   = contam_screen_hap1.output_fasta

        ## FCS outputs; hap2
        File fcs_hap2_int_clean_fa   = contam_screen_hap2.intermediate_clean_fa
        File fcs_hap2_contamFasta    = contam_screen_hap2.contamFasta
        File fcs_hap2_gx_report      = contam_screen_hap2.fcs_gx_report
        File fcs_hap2_tax_report     = contam_screen_hap2.fcs_taxonomy_report
        File fcs_hap2_adapter_report = contam_screen_hap2.adapter_Report
        File fcs_hap2_output_fasta   = contam_screen_hap2.output_fasta


        ## findMitosInHap1 outputs (mitos in assembly; not from mitoHiFi)
        File hap1_mt_contig_ids      = findMitosInHap1.mitoContigs
        File hap1_mt_blast           = findMitosInHap1.blastOutput
        File hap1_mt_parsed_blast    = findMitosInHap1.parsedBlast
        File hap1_mt_ncbi_blast      = findMitosInHap1.ncbiBlastOutput
        File hap1_mt_ncbi_filt_blast = findMitosInHap1.ncbiParsedBlast

        ## findMitosInHap2 outputs (mitos in assembly; not from mitoHiFi)
        File hap2_mt_contig_ids      = findMitosInHap2.mitoContigs
        File hap2_mt_blast           = findMitosInHap2.blastOutput
        File hap2_mt_parsed_blast    = findMitosInHap2.parsedBlast
        File hap2_mt_ncbi_blast      = findMitosInHap2.ncbiBlastOutput
        File hap2_mt_ncbi_filt_blast = findMitosInHap2.ncbiParsedBlast

        ## Assembly with mito contigs stripped out
        File hap1_mt_stripped_asm    = dropHap1Mito.FinalAssembly
        File hap2_mt_stripped_asm    = dropHap2Mito.FinalAssembly

        ## Assembly evaluation: hap1
        File hap1_T2Tcontigs         = evaluateHap1.T2Tcontigs
        File hap1_T2Tscaffolds       = evaluateHap1.T2Tscaffolds
        File hap1_unifiedAssembly    = evaluateHap1.unifiedAssembly
        File hap1_chr_names          = evaluateHap1.chromosome_names_from_CHM13
        File hap1_filteredFlanksBed  = evaluateHap1.filteredFlanksBed
        File hap1_bed_region         = evaluateHap1.bed_region
        File hap1_bed_SD             = evaluateHap1.bed_SD
        File hap1_bed_CENSAT         = evaluateHap1.bed_CENSAT
        File hap1_assemblyStatistics = evaluateHap1.assemblyStatistics

        ## Assembly evaluation: hap2
        File hap2_T2Tcontigs         = evaluateHap2.T2Tcontigs
        File hap2_T2Tscaffolds       = evaluateHap2.T2Tscaffolds
        File hap2_unifiedAssembly    = evaluateHap2.unifiedAssembly
        File hap2_chr_names          = evaluateHap2.chromosome_names_from_CHM13
        File hap2_filteredFlanksBed  = evaluateHap2.filteredFlanksBed
        File hap2_bed_region         = evaluateHap2.bed_region
        File hap2_bed_SD             = evaluateHap2.bed_SD
        File hap2_bed_CENSAT         = evaluateHap2.bed_CENSAT
        File hap2_assemblyStatistics = evaluateHap2.assemblyStatistics        
    }
}
