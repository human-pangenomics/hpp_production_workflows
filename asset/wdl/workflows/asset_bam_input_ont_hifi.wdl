version 1.0 

import "../tasks/asset.wdl" as asset_t
import "../tasks/bam2paf.wdl" as bam2paf_t
import "../tasks/bam_coverage.wdl" as bam_coverage_t
import "../tasks/tar.wdl" as tar_t

workflow assetTwoPlatforms {
    input {
        String sampleName
        String sampleSuffix
        Array[File] ontBamFiles
        Array[File] hifiBamFiles
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
    
    call bam_coverage_t.bamCoverage as ontBamCoverage {
        input:
            bamFiles = ontBamFiles,
            minMAPQ = minMAPQ
    }
 
    call bam_coverage_t.bamCoverage as hifiBamCoverage {
        input:
            bamFiles = hifiBamFiles,
            minMAPQ = minMAPQ
    }
    call asset_t.ast_pbTask as ontAssetTask{
        input:
            sampleName = "${sampleName}.${sampleSuffix}.ont",
            pafFiles = ontBam2Paf.pafFile,
            coverageMean = ontBamCoverage.coverageMean,
            coverageSD = ontBamCoverage.coverageSD,
            memSize = 32,
            threadCount = 8,
            diskSize = 256,
            preemptible = 2
    }
    call asset_t.ast_pbTask as hifiAssetTask{
        input:
            sampleName = "${sampleName}.${sampleSuffix}.hifi",
            pafFiles = hifiBam2Paf.pafFile,
            coverageMean = hifiBamCoverage.coverageMean,
            coverageSD = hifiBamCoverage.coverageSD,
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

