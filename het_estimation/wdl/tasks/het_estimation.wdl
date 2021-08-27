version 1.0 

import "../../../coverage/wdl/tasks/bam2paf.wdl" as bam2paf_t
import "../../../coverage/wdl/tasks/project_blocks.wdl" as project_t
import "../../../coverage/wdl/tasks/bedtools.wdl" as bedtools_t
import "../../../coverage/wdl/tasks/tar.wdl" as tar_t

workflow runHetEstimation{
    input {
        String sampleName
        File patAsm2RefBam
	File matAsm2RefBam
        File confidentRefBed
        File dipVcfGz
        File XYParBed
    }
    # Get heterozygousity statistics
    call hetEstimation {
        input:
            dipVcfGz = dipVcfGz,
            confidentRefBed = confidentRefBed
    }
    call bam2paf_t.bam2paf as patBam2Paf {
        input:
            bamFile = patAsm2RefBam,
            minMAPQ = 5,
            primaryOnly = "yes"
    }
    call bam2paf_t.bam2paf as matBam2Paf {
        input:
            bamFile = matAsm2RefBam,
            minMAPQ = 5,
            primaryOnly = "yes"
    }
    # Get three bed file pointing to the confident blocks in these three regions
    # 1. autosomal regions 
    # 2. PAR regions both X and Y
    # 3. non PAR regions just X
    call partitionBed {
        input:
            bed = confidentRefBed,
            XYParBed = XYParBed
    }
    # Project the three bed files to the assembly to find the 
    # assembly regions with confident alignments
    #pat
    call project_t.project as patProjectAutosome{
        input:
            blocksBed = partitionBed.autosomeBed,
            asm2refPaf = patBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "paternal.autosome",
            mode = "ref2asm"
    }
    call project_t.project as patProjectPar{
        input:
            blocksBed = partitionBed.parBed,
            asm2refPaf = patBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "paternal.par",
            mode = "ref2asm"
    }
    call project_t.project as patProjectNonParX{
        input:
            blocksBed = partitionBed.nonParChrXBed,
            asm2refPaf = patBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "paternal.non_par_chrX",
            mode = "ref2asm"
    }
    call tar_t.tarGz as patTarGz{
        input:
            tarGzName = "${sampleName}.paternal.projections",
            files = [patProjectAutosome.projectionBed, patProjectPar.projectionBed, patProjectNonParX.projectionBed]
    }
    # Project the three bed files to the assembly to find the
    # assembly regions with confident alignments
    #mat
    call project_t.project as matProjectAutosome{
        input:
            blocksBed = partitionBed.autosomeBed,
            asm2refPaf = matBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "maternal.autosome",
            mode = "ref2asm"
    }
    call project_t.project as matProjectPar{
        input:
            blocksBed = partitionBed.parBed,
            asm2refPaf = matBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "maternal.par",
            mode = "ref2asm"
    }

    call project_t.project as matProjectNonParX{
        input:
            blocksBed = partitionBed.parBed,
            asm2refPaf = matBam2Paf.pafFile,
            sampleName = sampleName,
            suffix = "maternal.par",
            mode = "ref2asm"
    }
    call tar_t.tarGz as matTarGz{
        input:
            tarGzName = "${sampleName}.maternal.projections",
            files = [matProjectAutosome.projectionBed, matProjectPar.projectionBed, matProjectNonParX.projectionBed]
    }
    output {
        File matProjectionsTarGz = matTarGz.fileTarGz
        File patProjectionsTarGz = patTarGz.fileTarGz
        File hetStatsText = hetEstimation.hetStatsText
    }
}

task partitionBed{
    input {
        File bed
        File XYParBed
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=16
        String dockerImage="quay.io/masri2019/hpp_base:latest"
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

        FILENAME=$(basename ~{bed})
        PREFIX=${FILENAME%.bed}

        mkdir outputs
        cat ~{bed} | awk '($1 != "chrY") && ($1 != "chrX")' > outputs/$PREFIX.autosome.bed
        cat ~{bed} | awk '$1 == "chrX"' > outputs/$PREFIX.chrX.bed
        bedtools intersect -a ~{bed} -b ~{XYParBed} > outputs/$PREFIX.par.bed
        bedtools subtract -a outputs/$PREFIX.chrX.bed -b ~{XYParBed} > outputs/$PREFIX.non_par_chrX.bed
        
    >>>
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File autosomeBed = glob("outputs/*.autosome.bed")[0]
        File parBed = glob("outputs/*.par.bed")[0]
        File nonParChrXBed = glob("outputs/*.non_par_chrX.bed")[0]
    }
}

task hetEstimation {
    input {
        File dipVcfGz
        File confidentRefBed
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=16
        String dockerImage="quay.io/masri2019/hpp_het_estiamtion:latest"
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
        
        FILENAME=$(basename ~{dipVcfGz})
        PREFIX=${FILENAME%.vcf.gz}
        
        bedtools intersect -header -a <(zcat ~{dipVcfGz}) -b ~{confidentRefBed} > $PREFIX.confident.vcf
        python3 ${HET_ESTIMATOR_PY} --vcf $PREFIX.confident.vcf > $PREFIX.het_stats.txt 
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File hetStatsText = glob("*.het_stats.txt")[0]
    }
}

