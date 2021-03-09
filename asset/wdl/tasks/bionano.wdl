version 1.0

import "tar.wdl" as tar_t

workflow runBionanoAlignment {
    input {
        String sampleName
        String sampleSuffix
        File sampleBnxGz
        File assemblyFastaGz
        Int preemptible=2
        String zones="us-west2-a"
    }
    call assembly2cmap {
        input:
            sampleBnxGz = sampleBnxGz,
            assemblyFastaGz = assemblyFastaGz,
            preemptible = preemptible,
            zones = zones
    }
    call alignBnx2Cmap {
        input:
            assemblyCmap = assembly2cmap.assemblyCmap,
            sampleBnxGz = sampleBnxGz,
            preemptible = preemptible,
            zones = zones
    }
    call tar_t.tarGz as bionanoTar{
        input:
            tarGzName = "${sampleName}.${sampleSuffix}.bionano.alignment",
            files = [alignBnx2Cmap.alignRefCmap, alignBnx2Cmap.alignQueryCmap, alignBnx2Cmap.alignXmap, assembly2cmap.assemblyKeyText]
    }
    output {
        File bionanoAlignmentTarGz = bionanoTar.fileTarGz
    }
}

task alignBnx2Cmap {
    input {
        File assemblyCmap
        File sampleBnxGz
        Int memSize=32
        Int threadCount=32
        Int diskSize=128
        String dockerImage="quay.io/masri2019/hpp_bionano:latest"
        Int preemptible=2
        String zones
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

        BNX_FILENAME=$(basename ~{sampleBnxGz})
        enzyme=`echo $BNX_FILENAME | cut -d_ -f3`

        gunzip -c ~{sampleBnxGz} > ${BNX_FILENAME%.gz}
        mkdir align
        python2 /home/apps/solve/Pipeline/1.0/align_bnx_to_cmap.py --prefix ${enzyme} --mol ${BNX_FILENAME%.gz} --ref ~{assemblyCmap} --ra /home/apps/solve/RefAligner/1.0 --nthreads ~{threadCount} --output align --optArgs /home/apps/solve/RefAligner/1.0/optArguments_nonhaplotype_DLE1_saphyr_human.xml --pipeline /home/apps/solve/Pipeline/1.0 
    >>>
    output {
        File alignRefCmap = glob("align/contigs/alignmolvref/merge/*_r.cmap")[0]
        File alignQueryCmap = glob("align/contigs/alignmolvref/merge/*_q.cmap")[0]
        File alignXmap = glob("align/contigs/alignmolvref/merge/*.xmap")[0]
    }
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
        continueOnReturnCode: true # Because of a warning for zipfile in python2.7 (zipfile.LargeZipFile: Filesize would require ZIP64 extensions)
    }

}

task assembly2cmap {
    input {
        File sampleBnxGz
        File assemblyFastaGz
        Int memSize=16
        Int threadCount=4
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_bionano:latest"
        Int preemptible=2
        String zones
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
        
        fn=`basename ~{sampleBnxGz}`
        fn_pref=`echo $fn | cut -d_ -f1`
        tech=`echo $fn | cut -d_ -f2`
        enzyme=`echo $fn | cut -d_ -f3`
        mkdir results
        FILENAME=$(basename "~{assemblyFastaGz}")
        PREFIX="${FILENAME%.fa.gz}" 
        gunzip -c ~{assemblyFastaGz} > ${PREFIX}.fa
        perl /home/apps/solve/fa2cmap_multi_color.pl -e ${enzyme:0:4} 1 -i ${PREFIX}.fa -o results
    >>>
    output {
        File assemblyCmap = glob("results/*.cmap")[0]
        File assemblyKeyText = glob("results/*_key.txt")[0]
        File fa2cmapSummaryText = glob("results/*_summary.txt")[0]
    }
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
        zones: zones
    }
 
}
