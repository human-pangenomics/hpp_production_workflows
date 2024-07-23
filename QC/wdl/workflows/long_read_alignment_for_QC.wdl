version 1.0

import "https://raw.githubusercontent.com/mobinasri/flagger/dev-0.3.0/wdls/workflows/long_read_aligner_scattered.wdl" as long_read_aligner
import "https://raw.githubusercontent.com/mobinasri/flagger/dev-0.3.0/wdls/tasks/other/misc.wdl" as misc_t

workflow long_read_alignment_for_qc_wf {

    input {
        File hap1Fasta # paternal
        File hap2Fasta # maternal
        String sampleName
        Array[File] hifi_reads
        Array[File] ont_reads

    }

    # Create a diploid assembly from hap1 & hap2
    call misc_t.createDipAsm {
        input:
            hap1AssemblyFasta   = hap1Fasta,
            hap2AssemblyFasta   = hap2Fasta,
            outputName          = "${sampleName}.dip.asm"
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

    output {

        ## minimap2 HiFi outputs
        File minimap2_hifi_bam            = align_hifi_minimap2.bamFile
        File minimap2_hifi_bai            = align_hifi_minimap2.baiFile
        File? minimap2_hifi_secphase_log  = align_hifi_minimap2.secphaseOutputLog

        ## minimap2 Ont outputs
        File minimap2_ont_bam             = align_ont_minimap2.bamFile
        File minimap2_ont_bai             = align_ont_minimap2.baiFile
        File? minimap2_ont_secphase_log   = align_ont_minimap2.secphaseOutputLog

   }
}