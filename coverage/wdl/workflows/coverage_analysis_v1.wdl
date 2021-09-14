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
        File filteredBedsTarGz = combineBeds.filteredBedsTarGz
    }
}

task combineBeds {
    input {
        File contigNamesText
        File genomeBedsTarGz
        File contigBedsTarGz
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
        
        mkdir genome_based contig_based
        tar --strip-components 1 -xvzf ~{genomeBedsTarGz} --directory genome_based
        tar --strip-components 1 -xvzf ~{contigBedsTarGz} --directory contig_based
                
        FILENAME=~{contigBedsTarGz}
        PREFIX=$(basename ${FILENAME%.*.*.tar.gz})
        
        mkdir genome_based_excluded
        grep -F -v -f ~{contigNamesText} genome_based/${PREFIX}.whole_genome_based.error.bed > genome_based_excluded/${PREFIX}.whole_genome_based.error.bed
        grep -F -v -f ~{contigNamesText} genome_based/${PREFIX}.whole_genome_based.duplicated.bed > genome_based_excluded/${PREFIX}.whole_genome_based.dup.bed
        grep -F -v -f ~{contigNamesText} genome_based/${PREFIX}.whole_genome_based.haploid.bed > genome_based_excluded/${PREFIX}.whole_genome_based.haploid.bed
        grep -F -v -f ~{contigNamesText} genome_based/${PREFIX}.whole_genome_based.collapsed.bed > genome_based_excluded/${PREFIX}.whole_genome_based.collapsed.bed

        mkdir combined
        cat genome_based_excluded/*.error.bed contig_based/*.error.bed | bedtools sort -i - > combined/${PREFIX}.combined.error.bed
        cat genome_based_excluded/*.dup.bed contig_based/*.duplicated.bed | bedtools sort -i - > combined/${PREFIX}.combined.duplicated.bed
        cat genome_based_excluded/*.haploid.bed contig_based/*.haploid.bed | bedtools sort -i - > combined/${PREFIX}.combined.haploid.bed
        cat genome_based_excluded/*.collapsed.bed contig_based/*.collapsed.bed | bedtools sort -i - > combined/${PREFIX}.combined.collapsed.bed

        mkdir filtered
        bedtools merge -d ~{mergeLength} -i combined/*.error.bed | awk '($3-$2) >= ~{minBlockLength}' > filtered/${PREFIX}.filtered.error.bed
        bedtools merge -d ~{mergeLength} -i combined/*.duplicated.bed | awk '($3-$2) >= ~{minBlockLength}' > filtered/${PREFIX}.filtered.duplicated.bed
        bedtools merge -d ~{mergeLength} -i combined/*.haploid.bed |  awk '($3-$2) >= ~{minBlockLength}' > filtered/${PREFIX}.filtered.haploid.bed
        bedtools merge -d ~{mergeLength} -i combined/*.collapsed.bed |  awk '($3-$2) >= ~{minBlockLength}' > filtered/${PREFIX}.filtered.collapsed.bed

        tar -cf ${PREFIX}.beds.combined.tar combined
        gzip ${PREFIX}.beds.combined.tar
        
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
        File combinedBedsTarGz = glob("*.beds.combined.tar.gz")[0]
        File filteredBedsTarGz = glob("*.beds.filtered.tar.gz")[0]
    }
}

