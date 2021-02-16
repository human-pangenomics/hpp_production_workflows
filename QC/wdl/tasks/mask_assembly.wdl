version 1.0

workflow runMaskAssembly {
    call maskAssembly
}

task maskAssembly {
    input {
        File assemblyFasta
        File adapterBed
        String maskedAssemblyName

        Int memSizeGB = 2
        Int threadCount = 2
        Int diskSizeGB = 64
        String dockerImage = "biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1"
    }

    command <<<

        ## Help w/ debugging...
        set -eux -o pipefail


        ## if bed file is empty (no adapters to mask): just copy to new file name
        if [[ ! -s adapterBed ]]; then
            cp ~{assemblyFasta} ~{maskedAssemblyName}.gz
        

        ## else, bed not empty (there are adapters to mask)
        else
            ## gunzip assembly if neccesary
            FILENAME=$(basename -- "~{assemblyFasta}")
            if [[ $FILENAME =~ \.gz$ ]]; then
                cp ~{assemblyFasta} .
                gunzip $FILENAME
            fi

            ## Get assembly filename (w/out .gz)
            ASM_ID=$(basename ~{assemblyFasta} | sed 's/.gz$//')

            bedtools maskfasta -fi ${ASM_ID} -bed ~{adapterBed} -fo ~{maskedAssemblyName}

            ## gzip assembly back up
            gzip ~{maskedAssemblyName}
        fi

    >>>

    output {
        File maskedAssembly = "~{maskedAssemblyName}.gz"

    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}