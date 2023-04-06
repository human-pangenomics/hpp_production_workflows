version 1.0 

workflow RunShortQC {
    input {
        File patFasta
        File matFasta
        File patYak
        File matYak
        File genesFasta
        File hs38Paf
        String childID
        # runtime configurations
        Int memSize
        Int threadCount
        String dockerImage = "quay.io/masri2019/qc-stats:latest"
        Int preemptible=0
        Int diskSize
        String zones="genomics.default-zones"
    }

    call shortQC{
        input:
            patFasta = patFasta ,
            matFasta = matFasta ,
            patYak = patYak ,
            matYak = matYak ,
            genesFasta = genesFasta ,
            hs38Paf = hs38Paf ,
            childID = childID ,
            # runtime configurations ,
            memSize = memSize ,
            threadCount = threadCount ,
            dockerImage = dockerImage ,
            preemptible = preemptible,
            diskSize = diskSize,
            zones = zones
    }

    output {
        File qcStatsText = shortQC.qcStatsText
    }
}

task shortQC {
    input{
        File patFasta
        File matFasta
        File patYak
        File matYak
        File childYak
        File genesFasta
        File hs38Paf
        String childID
        # runtime configurations    
        Int memSize
        Int threadCount
        String dockerImage
        Int preemptible
        Int diskSize
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

        # Computing length statistics
        k8 ${CAL_N50_PATH} ~{patFasta} > pat.len.stats.txt
        k8 ${CAL_N50_PATH} ~{matFasta} > mat.len.stats.txt

        # Aligning genes to assembly
        minimap2 -cx splice:hq -t ~{threadCount} ~{patFasta} ~{genesFasta} > pat.paf
        minimap2 -cx splice:hq -t ~{threadCount} ~{matFasta} ~{genesFasta} > mat.paf

        # hs38.ensembl.v99.cdna is already aligned to hg38 and its paf file should be given as an input
        # link to hs38.ensembl.v99.cdna : ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
        ln ~{hs38Paf} hs38.paf
        # Computing statistics for gene completeness
        paftools.js asmgene -a hs38.paf pat.paf > pat.gene.stats.txt
        paftools.js asmgene -a hs38.paf mat.paf > mat.gene.stats.txt

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{patFasta} > pat.error.stats.txt
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{matFasta} > mat.error.stats.txt

        # Computing QV stats
        yak qv -t 32 -p -K 3.2g -l 100k ~{childYak} <(zcat ~{patFasta}) > pat.yak_qv.txt
        yak qv -t 32 -p -K 3.2g -l 100k ~{childYak} <(zcat ~{matFasta}) > mat.yak_qv.txt
        
        # Merging all statistics in a single well-organized text file
        # ${QC_STATS_GENERATOR_PATH} is the path to a python script used for merging
        python3 ${QC_STATS_GENERATOR_PATH} --childID ~{childID} --patLenStats pat.len.stats.txt --matLenStats mat.len.stats.txt --patGeneStats pat.gene.stats.txt --matGeneStats mat.gene.stats.txt --patErrorStats pat.error.stats.txt --matErrorStats mat.error.stats.txt --patQVStats pat.yak_qv.txt --matQVStats mat.yak_qv.txt --output ~{childID}.qc.stats.txt

    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible: preemptible
        zones: zones
    }

    output {
        File qcStatsText = "~{childID}.qc.stats.txt"
    }
}
