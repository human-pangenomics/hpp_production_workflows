version 1.0 

import "../tasks/bwa.wdl" as bwa_t
import "../../../coverage/wdl/tasks/long_read_aligner.wdl" as longReadAligner_t

workflow fullAlignmentForAsset {
    input {
        String sampleName
        String sampleSuffix
        Array[File] ontReadFiles
        Array[File] hifiReadFiles
        Array[File] hicReadFiles_1
        Array[File] hicReadFiles_2
        File assembly
        String longReadAligner = "winnowmap"
        File? referenceFasta
        Int preemptible=2
        String zones="us-west2-a"
    }
    
    ## align HiC reads
    call bwa_t.bwaPairedAlignment as hicAlignment{
        input:
            sampleName = "${sampleName}.${sampleSuffix}.hic",
            readFiles_1 = hicReadFiles_1,
            readFiles_2 = hicReadFiles_2,
            assembly = assembly,
            referenceFasta = referenceFasta,
            preemptible = preemptible,
            zones = zones
    }

    ## align HiFi reads
    call longReadAligner_t.longReadAlignment as hifiAlignment {
        input:
            aligner = longReadAligner,
            preset = "map-pb",
            sampleName = "${sampleName}.${sampleSuffix}.hifi",
            readFiles = hifiReadFiles,
            assembly = assembly,
            referenceFasta = referenceFasta,
            preemptible = preemptible,
            zones = zones
    }

    ## align ONT reads
    call longReadAligner_t.longReadAlignment as ontAlignment {
        input:
            aligner = longReadAligner,
            preset = "map-ont",
            sampleName = "${sampleName}.${sampleSuffix}.ont",
            readFiles = ontReadFiles,
            assembly = assembly,
            referenceFasta = referenceFasta,
            preemptible = preemptible,
            zones = zones
    }

    output{
        ## HiFi outputs
        File hifiBam = hifiAlignment.sortedBamFile
        File hifiBai = hifiAlignment.baiFile
        ## ONT outputs
        File ontBam = ontAlignment.sortedBamFile
        File ontBai = ontAlignment.baiFile
        # HiC outputs
        File hicBam = hicAlignment.sortedBamFile
        File hicBai = hicAlignment.baiFile
    }
}

