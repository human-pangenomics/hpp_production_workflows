version 1.0 

workflow runStainedGlass{
    input {
        Array[File] veritymapTarGzArray
    }

    scatter (veritymapTarGz in veritymapTarGzArray) {

        # extract Fasta file
        call pullOutFasta{
            input:
                veritymapTarGz = veritymapTarGz
        }

        # Run StainedGlass for visualization
        call stainedGlass{
            input:
                cenFasta = pullOutFasta.fasta,
                cenFai = pullOutFasta.fai
        }
    }
    output{
        Array[File] stainedGlassTarGzArray = stainedGlass.outputTarGz
    }
}

task pullOutFasta {
    input {
        File veritymapTarGz
        # runtime configurations
        Int memSize=4
        Int threadCount=2
        Int diskSize=128
        String dockerImage="mobinasri/bio_base:latest"
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
        
        FILENAME=$(basename ~{veritymapTarGz})
        PREFIX=${FILENAME%%.tar.gz}

        mkdir output
        tar --strip-components 1 -xvzf ~{veritymapTarGz} --directory output
    >>>

    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File fasta = glob("output/*.fa")[0]
        File fai = glob("output/*.fai")[0]
    }
}

task stainedGlass {
    input {
       	File cenFasta
        File cenFai
        Int window = 2000
        Int mm_f = 10000
        # runtime configurations
        Int memSize=16
        Int threadCount=8
        Int diskSize=64
        String dockerImage="mobinasri/stained_glass:latest"
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

        FILENAME=$(basename ~{cenFasta})
        PREFIX=${FILENAME%%.fa*}


        WORK_DIR=$PWD
        cd /home/apps/StainedGlass
        ln -s ~{cenFasta} ${PREFIX}.fa
        ln -s ~{cenFai} ${PREFIX}.fa.fai
        printf "sample: ${PREFIX}\nfasta: ${PREFIX}.fa\nwindow: ~{window}\nnbatch: 1\nalnthreads: ~{threadCount}\nmm_f: ~{mm_f}\ntempdir: temp\n" > config/config.yaml
        conda config --set channel_priority strict
        conda run -n snakemake snakemake --use-conda --cores ~{threadCount} make_figures
          
        # Rename output folder
        mv results ${PREFIX}
        tar cvzf ${PREFIX}.tar.gz ${PREFIX}
        mv ${PREFIX}.tar.gz ${WORK_DIR}/
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File outputTarGz = glob("*.tar.gz")[0]
    }
}

