version 1.0

import "../tasks/cov2counts.wdl" as cov2counts_t
import "../tasks/cov2counts_contig_wise.wdl" as cov2counts_contig_wise_t
import "../tasks/fit_model.wdl" as fit_model_t
import "../tasks/fit_model_contig_wise.wdl" as fit_model_contig_wise_t
import "../tasks/find_blocks.wdl" as find_blocks_t
import "../tasks/find_blocks_contig_wise.wdl" as find_blocks_contig_wise_t
import "../tasks/pdf_generator.wdl" as pdf_generator_t
import "../tasks/bedtools.wdl" as bedtools_t
import "../tasks/fit_model_bed.wdl" as fit_model_bed_t

workflow runFlaggerPhase2{
    input {
        File hsatBedsTsv
        File coverageGz
        File highMapqCoverageGz
        File fai
        Float covFloat
        Boolean isDiploid
    }
    ## Each element in hsatBedsArray is an array itself;
    ## [BED URL, Coverage Factor, Suffix Name]
    Array[Array[String]] hsatBedsArray = read_tsv(hsatBedsTsv)
    scatter (hsatBed in hsatBedsArray){
        call bedtools_t.merge {
            input:
                bed = hsatBed[0],
                margin = 50000,
                outputPrefix = basename(hsatBed[0], ".bed")
        }
        call String2Float{
            input:
                str = hsatBed[1]
        }
        call fit_model_bed_t.runFitModelBed as hsatModels {
            input:
                bed = merge.mergedBed,
                suffix = hsatBed[2],
                coverageGz = coverageGz,
                covFloat = covFloat * String2Float.number
         }
    }
    call mergeHsatBeds {
        input:
            bedsTarGzArray = hsatModels.bedsTarGz
    }
    call cov2counts_t.cov2counts {
        input:
            coverageGz = coverageGz 
    }
    call fit_model_t.fitModel {
        input:
            counts = cov2counts.counts
    }
    call find_blocks_t.findBlocks {
        input:
            coverageGz = coverageGz,
            table = fitModel.probabilityTable
    }
    call cov2counts_contig_wise_t.cov2countsContigWise {
        input:
            coverageGz = coverageGz,
            fai = fai
    }
    call fit_model_contig_wise_t.fitModelContigWise {
        input:
            windowsText = cov2countsContigWise.windowsText,
            countsTarGz = cov2countsContigWise.contigCountsTarGz 
    }
    call find_blocks_contig_wise_t.findBlocksContigWise {
        input:
            contigCovsTarGz = cov2countsContigWise.contigCovsTarGz,
            contigProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz,
            windowsText = cov2countsContigWise.windowsText
    }
    call pdf_generator_t.pdfGenerator {
        input:
            contigProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz,
            genomeProbTable = fitModel.probabilityTable,
            isDiploid = isDiploid
    }
    call combineBeds as combineWindowBased{
        input:
            outputPrefix = "window_corrected",
            firstPrefix = "whole_genome",
            secondPrefix = "window_based",
            firstBedsTarGz = findBlocks.bedsTarGz,
            secondBedsTarGz = findBlocksContigWise.contigBedsTarGz
    }
    call combineBeds as combineHsatBased{
       input:
            outputPrefix = "hsat_corrected",
            firstPrefix = "window_corrected",
            secondPrefix = "hsat_based",
            firstBedsTarGz = combineWindowBased.combinedBedsTarGz,
            secondBedsTarGz = mergeHsatBeds.bedsTarGz
    }    
    call dupCorrectBeds {
        input:
            covGz = coverageGz,
            highMapqCovGz = highMapqCoverageGz,
            bedsTarGz = combineHsatBased.combinedBedsTarGz,
            prefix="hsat_corrected"
    }
    call filterBeds {
        input:
            fai = fai,
            dupCorrectedBedsTarGz = dupCorrectBeds.dupCorrectedBedsTarGz
    }
    output {
        File genomeCounts = cov2counts.counts
        File genomeProbTable = fitModel.probabilityTable
        File genomeBedsTarGz = findBlocks.bedsTarGz
        File windowCountsTarGz = cov2countsContigWise.contigCountsTarGz
        File windowCovsTarGz = cov2countsContigWise.contigCovsTarGz
        File windowProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz
        File windowBedsTarGz = findBlocksContigWise.contigBedsTarGz
        File pdf = pdfGenerator.pdf
        File combinedBedsTarGz = combineWindowBased.combinedBedsTarGz
        File dupCorrectedBedsTarGz = dupCorrectBeds.dupCorrectedBedsTarGz
        File filteredBedsTarGz = filterBeds.filteredBedsTarGz
        File hsatCorrectedBedsTarGz =  combineHsatBased.combinedBedsTarGz
    }
}

task String2Float {
    input {
        String str
    }
    command <<<
        echo ~{str} > file.txt
    >>>
    runtime {
        docker: "quay.io/masri2019/hpp_coverage:latest"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 1 SSD"
    }
    output {
        Float number = read_float("file.txt")
    }
}
task combineBeds {
    input {
        File firstBedsTarGz
        File secondBedsTarGz
        String firstPrefix
        String secondPrefix
        String outputPrefix = "combined"
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        
        mkdir first second
        tar --strip-components 1 -xvzf ~{firstBedsTarGz} --directory first
        tar --strip-components 1 -xvzf ~{secondBedsTarGz} --directory second
                
        FILENAME=~{firstBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})
        
        cat second/*.bed | sort -k1,1 -k2,2n | bedtools merge -i - > second_all.bed 
        mkdir first_minus_second ~{outputPrefix}
        for c in error duplicated haploid collapsed
        do
            bedtools subtract -sorted -a first/${PREFIX}.~{firstPrefix}.${c}.bed -b second_all.bed > first_minus_second/${PREFIX}.${c}.bed
            cat first_minus_second/*.${c}.bed second/*.${c}.bed | sort -k1,1 -k2,2n | bedtools merge -i - > ~{outputPrefix}/${PREFIX}.~{outputPrefix}.${c}.bed
        done

        tar -cf ${PREFIX}.beds.~{outputPrefix}.tar ~{outputPrefix}
        gzip ${PREFIX}.beds.~{outputPrefix}.tar

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File combinedBedsTarGz = glob("*.beds.${outputPrefix}.tar.gz")[0]
    }
}


task dupCorrectBeds {
    input {
        File covGz
        File highMapqCovGz
        File bedsTarGz
        String prefix
        Int minCov=4
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
 
        FILENAME=$(basename ~{highMapqCovGz})
        PREFIX=${FILENAME%.cov.gz}

        mkdir ~{prefix}
        tar --strip-components 1 -xvzf ~{bedsTarGz} --directory ~{prefix}

        zcat ~{highMapqCovGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= ~{minCov}) {print contig"\t"$1-1"\t"$2}}' | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - > high_mapq.bed


        zcat ~{covGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 < ~{minCov}) {print contig"\t"$1-1"\t"$2}}' | \
            sort -k1,1 -k2,2n | \
            bedtools merge -i - > extremely_low.bed

        mkdir dup_corrected

        cat high_mapq.bed extremely_low.bed | sort -k1,1 -k2,2n | bedtools merge -i - > exclude_dup.bed

        # do the correction
        bedtools subtract -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b exclude_dup.bed > dup_corrected/${PREFIX}.dup_corrected.duplicated.bed
        bedtools intersect -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b high_mapq.bed > dup_to_hap.bed
        ##bedtools intersect -sorted -a ~{prefix}/${PREFIX}.~{prefix}.duplicated.bed -b extremely_low.bed > dup_to_err.bed
        cat dup_to_hap.bed ~{prefix}/${PREFIX}.~{prefix}.haploid.bed | sort -k1,1 -k2,2n | bedtools merge -i - > dup_corrected/${PREFIX}.dup_corrected.haploid.bed
        cat extremely_low.bed ~{prefix}/${PREFIX}.~{prefix}.error.bed | sort -k1,1 -k2,2n | bedtools merge -i - > dup_corrected/${PREFIX}.dup_corrected.error.bed
        
        # just copy collapsed comp
        cp ~{prefix}/${PREFIX}.~{prefix}.collapsed.bed dup_corrected/${PREFIX}.dup_corrected.collapsed.bed

        tar -cf ${PREFIX}.beds.dup_corrected.tar dup_corrected
        gzip ${PREFIX}.beds.dup_corrected.tar

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File dupCorrectedBedsTarGz = glob("*.beds.dup_corrected.tar.gz")[0]
    }
}

task filterBeds {
    input {
        File fai
        File dupCorrectedBedsTarGz
        Int mergeLength=100
        Int minBlockLength=1000
        # runtime configurations
        Int memSize=8
        Int threadCount=4
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace
        
        mkdir dup_corrected
        tar --strip-components 1 -xvzf ~{dupCorrectedBedsTarGz} --directory dup_corrected

        FILENAME=~{dupCorrectedBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})

        mkdir initial filtered
        for c in error duplicated haploid collapsed
        do
            bedtools merge -d ~{mergeLength} -i dup_corrected/*.${c}.bed | awk '($3-$2) >= ~{minBlockLength}' > initial/${PREFIX}.filtered.${c}.bed
        done

        
        # Gather ambiguous overlaps
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.error.bed -b initial/${PREFIX}.filtered.haploid.bed > err_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.duplicated.bed -b initial/${PREFIX}.filtered.haploid.bed > dup_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.collapsed.bed -b initial/${PREFIX}.filtered.haploid.bed > col_hap.overlap.bed
        bedtools intersect -sorted -a initial/${PREFIX}.filtered.error.bed -b initial/${PREFIX}.filtered.duplicated.bed > err_dup.overlap.bed

        # Gather small blocks and assign as unknown
        cat initial/* | sort -k1,1 -k2,2n | bedtools merge -i - > initial/all.bed
        cat ~{fai} | awk '{print $1"\t"0"\t"$2}' | sort -k1,1 -k2,2n > asm.bed
        bedtools subtract -sorted -a asm.bed -b initial/all.bed > initial/${PREFIX}.filtered.unknown.bed
        # Add ambiguous overlaps as unknown
        cat *.overlap.bed initial/${PREFIX}.filtered.unknown.bed | sort -k1,1 -k2,2n | bedtools merge -i - > filtered/${PREFIX}.filtered.unknown.bed         

        # Subtract unknown regions
        for c in error duplicated haploid collapsed
        do
            bedtools subtract -sorted -a initial/${PREFIX}.filtered.${c}.bed -b filtered/${PREFIX}.filtered.unknown.bed > filtered/${PREFIX}.filtered.${c}.bed
        done

        tar -cf ${PREFIX}.beds.filtered.tar filtered
        gzip ${PREFIX}.beds.filtered.tar
        
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File filteredBedsTarGz = glob("*.beds.filtered.tar.gz")[0]
    }
}


task mergeHsatBeds {
    input {
        Array[File] bedsTarGzArray
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_coverage:latest"
        Int preemptible=2
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

         
        FILENAMES=(~{sep=" " bedsTarGzArray})
        FILENAME=${FILENAMES[0]}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})

        mkdir hsat_unmerged hsat_based
        for s in ~{sep=" " bedsTarGzArray}; do
            tar --strip-components 1 -xvzf $s --directory hsat_unmerged
        done
 
        for comp in error haploid duplicated collapsed; do
            cat hsat_unmerged/*.${comp}.bed | sort -k1,1 -k2,2n | bedtools merge -i - > hsat_based/$PREFIX.hsat_based.${comp}.bed
        done
         
        tar -cf ${PREFIX}.beds.hsat_based.tar hsat_based
        gzip ${PREFIX}.beds.hsat_based.tar
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File bedsTarGz = glob("*.beds.hsat_based.tar.gz")[0]
    }
}
