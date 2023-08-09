version 1.0

workflow RunFCS{
    input {
        File assembly
        
        File blast_div
        File GXI
        File GXS
        File manifest
        File metaJSON
        File seq_info
        File taxa

        Int threadCount
        Int preemptible = 1
        Int diskSizeGBGX  = 32
        Int diskSizeGBAdapter = 32

        String GxDB = basename(GXI, ".gxi")
        String asm_name=basename(sub(sub(assembly, "\\.gz$", ""), "\\.fasta$", ""))
    }

    call FCSGX {
        input:
            assembly=assembly,
            blast_div=blast_div,
            GXI=GXI,
            GXS=GXS,
            manifest=manifest,
            metaJSON=metaJSON,
            seq_info=seq_info,
            taxa=taxa,

            GxDB=GxDB,
            asm_name=asm_name,


            preemptible=preemptible,
            threadCount=threadCount,
            diskSizeGB=diskSizeGBGX
    }
    call FCS_adapter{
        input:
            cleanFasta = FCSGX.cleanFasta,
            asm_name=asm_name,


            preemptible=preemptible,
            threadCount=threadCount,
            diskSizeGB=diskSizeGBAdapter
    }

    output {
        File cleanFasta = FCSGX.cleanFasta
        File contamFasta = FCSGX.contamFasta
        File report = FCSGX.report
        
        File adapter_CleanedSequence = FCS_adapter.adapter_CleanedSequence
        File adapter_Report = FCS_adapter.adapter_Report
    }
    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

task FCSGX {
    input{
        File assembly
        File blast_div
        File GXI
        File GXS
        File manifest
        File metaJSON
        File seq_info
        File taxa
        

        String asm_name
        String GxDB

        Int memSizeGB = 500
        Int preemptible = 1
        Int diskSizeGB
        Int threadCount

    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        echo `df -h .`
        echo `du -h .`

        ln -s ~{blast_div}
        ln -s ~{GXI}
        ln -s ~{GXS}
        ln -s ~{manifest}
        ln -s ~{metaJSON}
        ln -s ~{seq_info}
        ln -s ~{taxa}

        ln -s ~{assembly}

        python3 /app/bin/run_gx --fasta ~{assembly} --gx-db ~{GxDB} --out-dir . --tax-id 9606
        zcat ~{assembly} | /app/bin/gx clean-genome --action-report ~{asm_name}.9606.fcs_gx_report.txt --output ~{asm_name}.clean.fasta --contam-fasta-out ~{asm_name}.contam.fasta 

        gzip ~{asm_name}.clean.fasta
        gzip ~{asm_name}.contam.fasta
    
    >>>

    output {
        File cleanFasta = "~{asm_name}.clean.fasta.gz"
        File contamFasta = "~{asm_name}.contam.fasta.gz"
        File report = "~{asm_name}.9606.fcs_gx_report.txt"
        
    }

    runtime {
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible : preemptible
        docker: 'ncbi/fcs-gx:0.4.0'
    }
}

task FCS_adapter {
    input {
        File cleanFasta

        String asm_name

        Int memSizeGB = 48
        Int preemptible = 1
        Int diskSizeGB
        Int threadCount
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # Run the adapter script 
        echo `df -h .`
        echo `du -h .`
        /app/fcs/bin/av_screen_x -o . --euk ~{cleanFasta}
        mv cleaned_sequences/* ~{asm_name}.adapterClean.fa # the output of FCS adapter is not actually gzipped

        rm -rf cleaned_sequences/

        gzip ~{asm_name}.adapterClean.fa

        
    >>>

    output {
        File adapter_CleanedSequence = "~{asm_name}.adapterClean.fa.gz"
        File adapter_Report = "fcs_adaptor_report.txt"
    }

    runtime {
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "ncbi/fcs-adaptor:0.4.0"
    }
}

