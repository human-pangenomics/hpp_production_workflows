import "https://raw.githubusercontent.com/mobinasri/flagger/dev-0.3.0/wdls/workflows/long_read_aligner_scattered.wdl" as long_read_aligner
import "https://raw.githubusercontent.com/mobinasri/flagger/dev-0.3.0/wdls/tasks/other/misc.wdl" as misc_t
import "https://raw.githubusercontent.com/kmiga/alphaAnnotation/main/cenSatAnnotation/centromereAnnotation.wdl" as censat_annotation
import ../tasks/compleasm.wdl as compleasm_wf
import ../tasks/nucfreq.wdl as nucfreq_wf

workflow assembly_qc_wf {

    input {
        File hap1Fasta # paternal
        File hap2Fasta # maternal
        String sampleName
        Boolean isMaleSample

        Array[File] hifi_reads
        Array[File] ont_reads

        File x_hap_compleasm_db
        File y_hap_compleasm_db
    }


    call compleasm_wf.compleasm as compleasm_hap1 {
        input:
            assembly     = hap1Fasta,
            lineage_tar  = if isMaleSample then y_hap_compleasm_db else x_hap_compleasm_db,
            lineage      = "primates"
    }

    call compleasm_wf.compleasm as compleasm_hap2 {
        input:
            assembly     = hap2Fasta,
            lineage_tar  = x_hap_compleasm_db,
            lineage      = "primates"
    }

    # Create a diploid assembly from hap1 & hap2
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFastaGz = hap1Fasta,
            hap2AssemblyFastaGz = hap2Fasta,
            outputName          = "${sampleName}.dip.asm"
    }

    # annotate censat regions
    call censat_annotation.centromereAnnotation as dip_censat_annotation {
        input:
            fasta = createDipAsm.diploidAssemblyFastaGz
    }


    call long_read_aligner.longReadAlignmentScattered as align_hifi_minimap2 {
        input:
            aligner               = "minimap2",
            preset                = "map-hifi",
            suffix                = "~{sampleName}_hifi_minimap2",
            alignerOptions        = "-I8g --eqx -Y -L --cs -p0.5 -y",
            readExtractionOptions = "-TMl,Mm",
            kmerSize              = 19,
            sampleName            = sampleName,
            assemblyFasta         = createDipAsm.diploidAssemblyFastaGz,
            readFiles             = hifi_reads,
            enableAddingMDTag     = true,
            enableRunningSecphase = true
    }

    call long_read_aligner.longReadAlignmentScattered as align_ont_minimap2 {
        input:
            aligner               = "minimap2",
            preset                = "map-ont",
            suffix                = "~{sampleName}_ont_minimap2",
            alignerOptions        = "-I8g --eqx -Y -L --cs -p0.5 -y",
            readExtractionOptions = "-TMl,Mm",
            kmerSize              = 17,
            sampleName            = sampleName,
            assemblyFasta         = createDipAsm.diploidAssemblyFastaGz,
            readFiles             = ont_reads,
            enableAddingMDTag     = true,
            enableRunningSecphase = true
    }

    call long_read_aligner.longReadAlignmentScattered as align_hifi_winnowmap {
        input:
            aligner               = "winnowmap",
            preset                = "map-pb",
            suffix                = "~{sampleName}_hifi_winnowmap",
            alignerOptions        = "-I8g --eqx -Y -L --cs -p0.5 -y",
            readExtractionOptions = "-TMl,Mm",
            kmerSize              = 15,
            sampleName            = sampleName,
            assemblyFasta         = createDipAsm.diploidAssemblyFastaGz,
            readFiles             = hifi_reads,
            enableAddingMDTag     = true,
            enableRunningSecphase = true
    }

    call long_read_aligner.longReadAlignmentScattered as align_ont_winnowmap {
        input:
            aligner               = "winnowmap",
            preset                = "map-ont",
            suffix                = "~{sampleName}_ont_winnowmap",
            alignerOptions        = "-I8g --eqx -Y -L --cs -p0.5 -y",
            readExtractionOptions = "-TMl,Mm",
            kmerSize              = 15,
            sampleName            = sampleName,
            assemblyFasta         = createDipAsm.diploidAssemblyFastaGz,
            readFiles             = ont_reads,
            enableAddingMDTag     = true,
            enableRunningSecphase = true
    }

    call nucfreq_wf as nucfreq_hifi_minimap2 {
        input:
            input_bam             = align_hifi_minimap2.bamFile
            input_bam_bai         = align_hifi_minimap2.baiFile
            regions_bed           = dip_censat_annotation.centromeres
            assembly_fasta        = createDipAsm.diploidAssemblyFastaGz
    }

    call nucfreq_wf as nucfreq_ont_minimap2 {
        input:
            input_bam             = align_ont_minimap2.bamFile
            input_bam_bai         = align_ont_minimap2.baiFile
            regions_bed           = dip_censat_annotation.centromeres
            assembly_fasta        = createDipAsm.diploidAssemblyFastaGz
    }

    call nucfreq_wf as nucfreq_hifi_winnowmap {
        input:
            input_bam             = align_hifi_winnowmap.bamFile
            input_bam_bai         = align_hifi_winnowmap.baiFile
            regions_bed           = dip_censat_annotation.centromeres
            assembly_fasta        = createDipAsm.diploidAssemblyFastaGz
    }

    call nucfreq_wf as nucfreq_ont_winnowmap {
        input:
            input_bam             = align_ont_winnowmap.bamFile
            input_bam_bai         = align_ont_winnowmap.baiFile
            regions_bed           = dip_censat_annotation.centromeres
            assembly_fasta        = createDipAsm.diploidAssemblyFastaGz
    }


    output {

        ## hap1 compleasm
        File compleasm_summary_hap1      = compleasm_hap1.summary
        File compleasm_full_table_hap1   = compleasm_hap1.fullTable
        File compleasm_tar_hap1          = compleasm_hap1.outputTar

        ## hap2 compleasm
        File compleasm_summary_hap2      = compleasm_hap2.summary
        File compleasm_full_table_hap2   = compleasm_hap2.fullTable
        File compleasm_tar_hap2          = compleasm_hap2.outputTar

        
        ## censat outputs
        File as_hor_sf_bed                = dip_censat_annotation.as_hor_sf_bed
        File as_strand_bed                = dip_censat_annotation.as_strand_bed
        File as_hor_bed                   = dip_censat_annotation.as_hor_bed
        File as_sf_bed                    = dip_censat_annotation.as_sf_bed
        File cenSatAnnotations            = dip_censat_annotation.cenSatAnnotations
        File centromeres                  = dip_censat_annotation.centromeres


        ## minimap2 HiFi outputs
        File minimap2_hifi_bam            = align_hifi_minimap2.bamFile
        File minimap2_hifi_bai            = align_hifi_minimap2.baiFile
        File? minimap2_hifi_secphase_log  = align_hifi_minimap2.secphaseOutputLog

        ## minimap2 Ont outputs
        File minimap2_ont_bam             = align_ont_minimap2.bamFile
        File minimap2_ont_bai             = align_ont_minimap2.baiFile
        File? minimap2_ont_secphase_log   = align_ont_minimap2.secphaseOutputLog

        ## winnowmap HiFi outputs
        File winnowmap_hifi_bam            = align_hifi_winnowmap.bamFile
        File winnowmap_hifi_bai            = align_hifi_winnowmap.baiFile
        File? winnowmap_hifi_secphase_log  = align_hifi_winnowmap.secphaseOutputLog

        ## winnowmap Ont outputs
        File winnowmap_ont_bam             = align_ont_winnowmap.bamFile
        File winnowmap_ont_bai             = align_ont_winnowmap.baiFile
        File? winnowmap_ont_secphase_log   = align_ont_winnowmap.secphaseOutputLog


        ## nucfreq minimap2 hifi
        File nucplot_images_mm2_hifi       = nucfreq_hifi_minimap2.nucplot_image_tar
        File nucfreq_all_bed_mm2_hifi      = nucfreq_hifi_minimap2.nucfreq_all_bed
        File error_clusters_bed_mm2_hifi   = nucfreq_hifi_minimap2.error_clusters_bed      
        File first_allele_bigwig_mm2_hifi  = nucfreq_hifi_minimap2.first_allele_bigwig
        File second_allele_bigwig_mm2_hifi = nucfreq_hifi_minimap2.second_allele_bigwig

        ## nucfreq minimap2 ont
        File nucplot_images_mm2_ont        = nucfreq_ont_minimap2.nucplot_image_tar
        File nucfreq_all_bed_mm2_ont       = nucfreq_ont_minimap2.nucfreq_all_bed
        File error_clusters_bed_mm2_ont    = nucfreq_ont_minimap2.error_clusters_bed      
        File first_allele_bigwig_mm2_ont   = nucfreq_ont_minimap2.first_allele_bigwig
        File second_allele_bigwig_mm2_ont  = nucfreq_ont_minimap2.second_allele_bigwig

        ## nucfreq winnowmap hifi
        File nucplot_images_wm_hifi        = nucfreq_hifi_winnowmap.nucplot_image_tar
        File nucfreq_all_bed_wm_hifi       = nucfreq_hifi_winnowmap.nucfreq_all_bed
        File error_clusters_bed_wm_hifi    = nucfreq_hifi_winnowmap.error_clusters_bed      
        File first_allele_bigwig_wm_hifi   = nucfreq_hifi_winnowmap.first_allele_bigwig
        File second_allele_bigwig_wm_hifi  = nucfreq_hifi_winnowmap.second_allele_bigwig

        ## nucfreq winnowmap ont
        File nucplot_images_wm_ont         = nucfreq_ont_winnowmap.nucplot_image_tar
        File nucfreq_all_bed_wm_ont        = nucfreq_ont_winnowmap.nucfreq_all_bed
        File error_clusters_bed_wm_ont     = nucfreq_ont_winnowmap.error_clusters_bed      
        File first_allele_bigwig_wm_ont    = nucfreq_ont_winnowmap.first_allele_bigwig
        File second_allele_bigwig_wm_ont   = nucfreq_ont_winnowmap.second_allele_bigwig
    }
}