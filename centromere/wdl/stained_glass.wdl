version 1.0 

workflow runStainedGlass{
    input {
        Array[File] horAssemblyFastaArray
    }

    scatter (fasta in horAssemblyFastaArray) {
        # Run StainedGlass for visualization
        call stainedGlass{
            input:
                cenFasta = fasta
        }
    }
    output{
        Array[File] stainedGlassTarGzArray = stainedGlass.outputTarGz
    }
}


task stainedGlass {
    input {
       	File cenFasta
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


        # copy StainedGlass dir to the working dir
        cp -r /home/apps/StainedGlass .
        cd StainedGlass
        cp ~{cenFasta} ${PREFIX}.fa
        # index fasta
        samtools faidx ${PREFIX}.fa
        # configure inputs for running StainedGlass
        printf "sample: ${PREFIX}\nfasta: ${PREFIX}.fa\nwindow: ~{window}\nnbatch: 1\nalnthreads: ~{threadCount}\nmm_f: ~{mm_f}\ntempdir: temp\n" > config/config.yaml
        conda config --set channel_priority strict
        # run StainedGlass
        conda run -n snakemake snakemake --use-conda --cores ~{threadCount} make_figures
          
        # move results out of StainedGlass dir
        mv results ../${PREFIX}
        # exit the StainedGlass dir
        cd ../
        # compress results
        tar cvzf ${PREFIX}.tar.gz ${PREFIX}
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

