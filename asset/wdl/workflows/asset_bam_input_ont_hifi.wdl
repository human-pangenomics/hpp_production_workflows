version 1.0

import "../tasks/asset.wdl" as asset_t
import "../../../coverage/wdl/tasks/bam2paf.wdl" as bam2paf_t
import "../../../coverage/wdl/tasks/bam_coverage.wdl" as bam_coverage_t
import "../tasks/tar.wdl" as tar_t

workflow assetTwoPlatforms {
    input {
        String sampleName
        String sampleSuffix
        Array[File] ontBamFiles
        Array[File] hifiBamFiles
        Float ontCoverageMean
        Float ontCoverageSd
        Float hifiCoverageMean
        Float hifiCoverageSd
        Int minMAPQ = 21
    }
    
    scatter(bamFile in ontBamFiles) {
        call bam2paf_t.bam2paf as ontBam2Paf {
            input:
                bamFile = bamFile,
                minMAPQ = minMAPQ
        }
    }

    scatter(bamFile in hifiBamFiles) {
        call bam2paf_t.bam2paf as hifiBam2Paf {
            input: 
                bamFile = bamFile,
                minMAPQ = minMAPQ
        }
    }
    
    call asset_t.ast_pbTask as ontAssetTask{
        input:
            sampleName = "${sampleName}.${sampleSuffix}.ont",
            pafFiles = ontBam2Paf.pafFile,
            coverageMean = ontCoverageMean,
            coverageSD = ontCoverageSd,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = 2
    }
    call asset_t.ast_pbTask as hifiAssetTask{
        input:
            sampleName = "${sampleName}.${sampleSuffix}.hifi",
            pafFiles = hifiBam2Paf.pafFile,
            coverageMean = hifiCoverageMean,
            coverageSD = hifiCoverageSd,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = 2
    }

    call tar_t.tarGz as assetTar{
        input:
            tarGzName = "${sampleName}.${sampleSuffix}.asset",
            files = [hifiAssetTask.supportBed, ontAssetTask.supportBed]
    }
    output {
        File assetOutputsTarGz = assetTar.fileTarGz
    }
}

