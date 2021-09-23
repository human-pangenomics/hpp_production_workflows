version 1.0

import "../tasks/cov2counts.wdl" as cov2counts_t
import "../tasks/cov2counts_contig_wise.wdl" as cov2counts_contig_wise_t
import "../tasks/fit_model.wdl" as fit_model_t
import "../tasks/fit_model_contig_wise.wdl" as fit_model_contig_wise_t
import "../tasks/find_blocks.wdl" as find_blocks_t
import "../tasks/find_blocks_contig_wise.wdl" as find_blocks_contig_wise_t
import "../tasks/pdf_generator.wdl" as pdf_generator_t

workflow runCoverageAnalysisV1{
    input {
        File coverageGz
        File highMapqCoverageGz
        File fai
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
            coverageGz = coverageGz
    }
    call fit_model_contig_wise_t.fitModelContigWise {
        input:
            fai = fai,
            countsTarGz = cov2countsContigWise.contigCountsTarGz 
    }
    call find_blocks_contig_wise_t.findBlocksContigWise {
        input:
            contigCovsTarGz = cov2countsContigWise.contigCovsTarGz,
            contigProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz,
            contigNames = fitModelContigWise.longContigNamesText 
    }
    call pdf_generator_t.pdfGenerator {
        input:
            contigProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz,
            genomeProbTable = fitModel.probabilityTable
    }
    call combineBeds {
        input:
            contigNamesText = fitModelContigWise.longContigNamesText,
            genomeBedsTarGz = findBlocks.bedsTarGz,
            contigBedsTarGz = findBlocksContigWise.contigBedsTarGz
    }
    
    call dupCorrectBeds {
        input:
            highMapqCovGz = highMapqCoverageGz,
            combinedBedsTarGz = combineBeds.combinedBedsTarGz
    }
    call filterBeds {
        input:
            dupCorrectedBedsTarGz = dupCorrectBeds.dupCorrectedBedsTarGz
    }
    output {
        File genomeCounts = cov2counts.counts
        File genomeProbTable = fitModel.probabilityTable
        File genomeBedsTarGz = findBlocks.bedsTarGz
        File contigCountsTarGz = cov2countsContigWise.contigCountsTarGz
        File contigCovsTarGz = cov2countsContigWise.contigCovsTarGz
        File contigProbTablesTarGz = fitModelContigWise.contigProbTablesTarGz
        File contigBedsTarGz = findBlocksContigWise.contigBedsTarGz
        File pdf = pdfGenerator.pdf
        File combinedBedsTarGz = combineBeds.combinedBedsTarGz
        File dupCorrectedBedsTarGz = dupCorrectBeds.dupCorrectedBedsTarGz
        File filteredBedsTarGz = filterBeds.filteredBedsTarGz
    }
}

task combineBeds {
    input {
        File contigNamesText
        File genomeBedsTarGz
        File contigBedsTarGz
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
        
        mkdir genome_based contig_based
        tar --strip-components 1 -xvzf ~{genomeBedsTarGz} --directory genome_based
        tar --strip-components 1 -xvzf ~{contigBedsTarGz} --directory contig_based
                
        FILENAME=~{contigBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})
        
        mkdir genome_based_excluded combined
        for c in error duplicated haploid collapsed
        do
            if [ -s "genome_based/${PREFIX}.whole_genome_based.${c}.bed" ]
            then
                grep -F -v -f ~{contigNamesText} genome_based/${PREFIX}.whole_genome_based.${c}.bed > genome_based_excluded/${PREFIX}.whole_genome_based.${c}.bed
            else
                echo "" > genome_based_excluded/${PREFIX}.whole_genome_based.${c}.bed
            fi
            cat genome_based_excluded/*.${c}.bed contig_based/*.${c}.bed | bedtools sort -i - | bedtools merge -i - > combined/${PREFIX}.combined.${c}.bed
        done

        tar -cf ${PREFIX}.beds.combined.tar combined
        gzip ${PREFIX}.beds.combined.tar

    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File combinedBedsTarGz = glob("*.beds.combined.tar.gz")[0]
    }
}


task dupCorrectBeds {
    input {
        File highMapqCovGz
        File combinedBedsTarGz
        Int minCov=5
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

        combined
        tar --strip-components 1 -xvzf ~{combinedBedsTarGz} --directory combined

        zcat ~{highMapqCovGz} | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= ~{minCov}) {print contig"\t"$1-1"\t"$2}}' | \
            bedtools sort -i - | bedtools merge -i - > high_mapq.bed

        mkdir dup_corrected

        # do the correction
        bedtools subtract -a combined/${PREFIX}.duplicated.bed -b high_mapq.bed > dup_corrected/${PREFIX}.dup_corrected.duplicated.bed
        bedtools intersect -a combined/${PREFIX}.duplicated.bed -b high_mapq.bed > dup_to_hap.bed
        cat dup_to_hap.bed combined/${PREFIX}.haploid.bed | bedtools sort -i - | bedtools merge -i - > dup_corrected/${PREFIX}.dup_corrected.haploid.bed
        
        # just copy error and collapsed comps
        cp combined/${PREFIX}.combined.error.bed dup_corrected/${PREFIX}.dup_corrected.error.bed
        cp combined/${PREFIX}.combined.collapsed.bed dup_corrected/${PREFIX}.dup_corrected.collapsed.bed

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
        File dupCorrectedBedsTarGz
        Int mergeLength=100
        Int minBlockLength=1000
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
        
        mkdir dup_corrected
        tar --strip-components 1 -xvzf ~{dupCorrectedBedsTarGz} --directory dup_corrected

        FILENAME=~{dupCorrectedBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})

        mkdir filtered
        for c in error duplicated haploid collapsed
        do
            bedtools merge -d ~{mergeLength} -i dup_corrected/*.${c}.bed | awk '($3-$2) >= ~{minBlockLength}' > filtered/${PREFIX}.filtered.${c}.bed
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

