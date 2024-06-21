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
        lineage_tar: "BUSCO lineage downloaded with compleasm download command. Must match the lineage provided. Cannot use vanilla BUSCO lineages."
        otherArgs: "(default is empty string) Arguments to be passed to compleasm run command such as '--min_complete 0.9'"
        compleasm_script: "(optional) If passed use this script to run compleasm. Useful for updates that haven't been pushed to the official Docker image."
    }
  
    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail
        
        ## extract lineage_tar to the expected directory
        mkdir mb_download
        cd mb_download

        tar -zxvf ~{lineage_tar}
        touch ~{lineage}.done

        ## Check that extracted lineage matches the lineage specified 
        if [[ ! -d "~{lineage}_odb10" && ! -d "~{lineage}" ]]; then
            echo "lineage tar must create a folder with same name as lineage provided"
            exit 1
        fi

        cd ..

        ## get assembly prefix, remove fa, fa.gz, fasta, fasta.gz suffixes
        FILEPREFIX=$(basename ~{assembly} | sed 's/\(.*\)\(\.fa\|\.fa\.gz\|\.fasta\|\.fasta\.gz\)$/\1/')


        ## Actual run: miniprot --> hmm filter --> summarize
        ## This command includes the lineage folder due to how long the download takes using
        ## the compleasm download/run commands
        if [ "~{compleasm_script##*.}" = "py" ]; then
            python "~{compleasm_script}" run \
                -a ~{assembly} \
                -o ${FILEPREFIX}_compleasm \
                -t ~{threadCount} \
                -l ~{lineage} \
                -L ./mb_download \
                ~{otherArgs}
        else
            compleasm run \
                -a ~{assembly} \
                -o ${FILEPREFIX}_compleasm \
                -t ~{threadCount} \
                -l ~{lineage} \
                -L ./mb_download \
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
