version 1.0

# This is a task level wdl workflow to run the program COMPLEASM

workflow runCompleasm {
    call compleasm
  
    output {
        File compleasmSummary   = compleasm.summary
        File compleasmFullTable = compleasm.fullTable
        File compleasmOutputTar = compleasm.outputTar
    }

    meta {
        description: "Calls compleasm run to assess BUSCO gene completeness of an assembly. See [compleasm documentation](https://github.com/huangnengCSU/compleasm)"
    }
}

task compleasm {
    input{
        File assembly
        File lineage_tar
        String lineage     = "primates"
        String otherArgs   = ""

        File? compleasm_script

        Int threadCount    = 24    
        Int memSizeGB      = 48
        Int diskSizeGB     = 64    
        String dockerImage = "huangnengcsu/compleasm@sha256:34643e3ebbc9a2dde9d648434e1b79cc576b38967cad86341622de5701da8ede" # v0.2.6
    }

    parameter_meta {
        assembly: "Assembly to analyze for gene completeness in BUSCO genes."
        lineage: "(default is primates) Busco lineage name. For human use primates. This parameter must match the name of the folder that is extracted from the lineage_tar."
        lineage_tar: "BUSCO lineage download folder -- output of compleasm download command or prior compleasm run. Must match the lineage provided. Cannot use vanilla BUSCO lineages."
        otherArgs: "(default is empty string) Arguments to be passed to compleasm run command such as '--min_complete 0.9'"
        compleasm_script: "(optional) If passed use this script to run compleasm. Useful for updates that haven't been pushed to the official Docker image."
    }
  
    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        
        ## Find name of the lineage download folder. It is mb_download by convention
        ## but could be named something else. This is the folder which lineages are 
        ## downloaded into, so a primate lineage will also contain a eukaryote lineage, for example.
        LINEAGE_DOWNLOAD_DIR=$(tar -tzf ~{lineage_tar} | head -1 | cut -d/ -f1)

        ## Extract lineage_tar. This should be the entire download directory from a prior compleasm run
        ## or compleasm download command.
        tar -zvf ~{lineage_tar}


        ## Check that extracted lineage matches the lineage specified 
        if [[ ! -d "${LINEAGE_DOWNLOAD_DIR}/~{lineage}_odb10" && ! -d "${LINEAGE_DOWNLOAD_DIR}/~{lineage}" ]]; then
            echo "lineage tar must create a folder with same name as lineage provided"
            exit 1
        fi

        ## get assembly prefix, remove fa, fa.gz, fasta, fasta.gz suffixes
        FILEPREFIX=$(basename ~{assembly} | sed 's/\(.*\)\(\.fa\|\.fa\.gz\|\.fasta\|\.fasta\.gz\)$/\1/')


        ## Actual run: miniprot --> hmm filter --> summarize
        ## This command includes the lineage folder due to how long the download takes using
        ## the compleasm download/run commands
        
        ## put optional updated script into bash variable so we can check if it exists
        UPDATED_SCRIPT="~{compleasm_script}"

        if [ "${UPDATED_SCRIPT##*.}" = "py" ]; then
            python "~{compleasm_script}" run \
                -a ~{assembly} \
                -o ${FILEPREFIX}_compleasm \
                -t ~{threadCount} \
                -l ~{lineage} \
                -L "$LINEAGE_DOWNLOAD_DIR" \
                ~{otherArgs}
        else
            compleasm run \
                -a ~{assembly} \
                -o ${FILEPREFIX}_compleasm \
                -t ~{threadCount} \
                -l ~{lineage} \
                -L "$LINEAGE_DOWNLOAD_DIR" \
                ~{otherArgs}
        fi

        ## copy and rename summary file for delocalization (keep original in place)
        cp ${FILEPREFIX}_compleasm/summary.txt ${FILEPREFIX}_compleasm.summary.txt

        ## copy and rename full table for delocalization (keep original in place)
        cp ${FILEPREFIX}_compleasm/*/full_table.tsv ${FILEPREFIX}_compleasm.full_table.txt

        ## tar.gz results (including hmmer_output and miniprot_output.gff)
        tar -czf ${FILEPREFIX}_compleasm.tar.gz ${FILEPREFIX}_compleasm

  >>>

  output {
    File summary   = glob("*_compleasm.summary.txt")[0]
    File fullTable = glob("*_compleasm.full_table.txt")[0]
    File outputTar = glob("*_compleasm.tar.gz")[0]
  }

  runtime {
    memory: memSizeGB + " GB"
    cpu: threadCount
    disks: "local-disk " + diskSizeGB + " SSD"
    docker: dockerImage
    preemptible: 1
  }
}
