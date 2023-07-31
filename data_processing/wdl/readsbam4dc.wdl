version 1.0

workflow filter_reads_bam_wf {
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Download reads.bam files from s3 bucket and Q9 filter then demultiplex."
    }

    parameter_meta {
        s3_uris:       "Files should be named s3://blah/blah.blah"
        barcode:       "Barcode that the sample you are targeting uses -- such as bc1019."
        barcodes_file: "See: https://www.pacb.com/wp-content/uploads/Sequel_96_barcodes_v1.zip"
    }

    input {
        String s3_uri
        String barcode
        File barcodes_file
    }


    call filter_reads_bam {
        input:
            s3_uri        = s3_uri,
            barcode       = barcode,
            barcodes_file = barcodes_file
    }


   output {
        File q9_demux_bam  = filter_reads_bam.q9_demux_bam
        File lima_counts   = filter_reads_bam.lima_counts
        File lima_summary  = filter_reads_bam.lima_summary
    }  
}

task filter_reads_bam {
    input {
        String s3_uri
        String barcode
        File barcodes_file

        Int memSizeGB   = 16
        Int threadCount = 4
        Int diskSize    = 500
        Int preempts    = 3
    }


    command <<<

        set -eux -o pipefail

        mkdir download

        ## download file from S3 bucket
        ./s5cmd --no-sign-request cp  \
            ~{s3_uri} \
            download/

        ## pull reads over q9, remove kinetics tags (that make the file big)
        samtools view \
            --bam \
            --with-header \
            --remove-tag 'fi,ri,fp,rp' \
            --expr '[rq]>=0.88' \
            download/*.bam \
            q9.bam   

        input_bam_prefix=$(basename ~{s3_uri} .reads.bam)

        ## demultiplex
        lima \
            q9.bam \
            ~{barcodes_file} \
            ${input_bam_prefix}.out.bam \
            --dump-removed \
            --split-named \
            --min-ref-span 0.75 \
            --min-score 70 \
            --same \
            --single-side  

    >>>

    output {
        File q9_demux_bam  = glob("*~{barcode}.bam")[0]
        File lima_counts   = glob("*.out.lima.counts")[0]
        File lima_summary  = glob("*.out.lima.summary")[0]
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "quay.io/biocontainers/lima:2.7.1--h9ee0642_0"
        preemptible: preempts
    }
}
