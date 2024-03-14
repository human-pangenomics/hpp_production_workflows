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

        Int threadCount    = 24    
        Int memSizeGB      = 48
        Int diskSizeGB     = 64    
        String dockerImage = "huangnengcsu/compleasm@sha256:795bcdceaa25f15c14b11f2b8c5034423e5efc3310003746beb2d3ae09194d7d" # v0.2.5
    }

    parameter_meta {
        assembly: "Assembly to analyze for gene completeness in BUSCO genes."
        lineage: "(default is primates) Busco lineage name. For human use primates."
        lineage_tar: "BUSCO lineage downloaded with compleasm download command. Must match the lineage provided. Cannot use vanilla BUSCO lineages."        
        otherArgs: "(default is empty string) Arguments to be passed to compleasm run command such as '--min_complete 0.9'"
    }
  
    command <<<

        # exit when a command fails, fail with unset variables, print commands before execution
        set -eux -o pipefail

        ## extract lineage_tar
        mkdir lineage_folder
        tar -xzf ~{lineage_tar} -C lineage_folder

        ## get assembly prefix, remove fa, fa.gz, fasta, fasta.gz suffixes
        FILEPREFIX=$(basename ~{assembly} | sed 's/\(.*\)\(\.fa\|\.fa\.gz\|\.fasta\|\.fasta\.gz\)$/\1/')


        ## Actual run: miniprot --> hmm filter --> summarize
        ## This command includes the lineage folder due to how long the download takes using
        ## the compleasm download/run commands
        compleasm run \
          -a ~{assembly} \
          -o ${FILEPREFIX}_compleasm \
          -t ~{threadCount} \
          -l ~{lineage} \
          -L ./lineage_folder \
          ~{otherArgs}


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
