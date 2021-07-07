version 1.0 

import "../tasks/asset.wdl" as asset_t
import "../../../coverage/wdl/tasks/bam2paf.wdl" as bam2paf_t
import "../../../coverage/wdl/tasks/bam_coverage.wdl" as bam_coverage_t
import "../tasks/tar.wdl" as tar_t

workflow assetFourPlatforms {
    input {
        String sampleName
        String sampleSuffix
        Array[File] ontBamFiles
        Array[File] hifiBamFiles
        Array[File] hicBamFiles
        File bionanoAlignmentTarGz
        File assembly
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
            sampleName = "${sampleName}.${sampleSuffix}.ont",
            bamFiles = ontBamFiles,
            minMAPQ = minMAPQ
    }
 
    call bam_coverage_t.bamCoverage as hifiBamCoverage {
        input:
            sampleName = "${sampleName}.${sampleSuffix}.hifi",
            bamFiles = hifiBamFiles,
            minMAPQ = minMAPQ
    }
    call asset_t.asset as asset {
        input:
            sampleName = "${sampleName}.${sampleSuffix}",
            ontPafFiles = ontBam2Paf.pafFile,
            hifiPafFiles = hifiBam2Paf.pafFile,
            hicBamFiles = hicBamFiles,
            bionanoAlignmentTarGz = bionanoAlignmentTarGz,
            assembly = assembly,
            ontCoverageMean = ontBamCoverage.coverageMean, 
            ontCoverageSD = ontBamCoverage.coverageSD,
            hifiCoverageMean = hifiBamCoverage.coverageMean,
            hifiCoverageSD = hifiBamCoverage.coverageSD
    }
    call tar_t.tarGz as assetTar{
        input:
            tarGzName = "${sampleName}.${sampleSuffix}.asset",
            files = [asset.gapsBed, asset.hicBed, asset.hifiBed, asset.ontBed, asset.bionanoBed]
    }
    call tar_t.tarGz as coverageFilesTar{
        input:
            tarGzName = "${sampleName}.${sampleSuffix}.coverage",
            files = [asset.hicCoverageWig, asset.hifiCoverageWig, asset.ontCoverageWig, asset.bionanoCoverageWig, asset.hicCoverage, hifiBamCoverage.coverageGz, ontBamCoverage.coverageGz]
    }
    output {
        File assetOutputsTarGz = assetTar.fileTarGz
        File coverageFilesTarGz = coverageFilesTar.fileTarGz
    }
}

