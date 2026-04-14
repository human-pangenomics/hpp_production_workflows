version 1.0


task ntsm_count {
    input {
        Array[File] input_reads
        String sample_id
        String read_type

        File? cram_reference

        Int max_coverage = 10
        Int memSizeGB    = 4
        Int threadCount  = 4
        Int addldisk     = 10
        Int preempts     = 2
    }

    parameter_meta {
        input_reads:  "Files must be in fastq, fastq.gz (preferred), BAM, or CRAM format. CRAM files must be passed with cram_reference"
        max_coverage: "Stop processing once average site coverage reaches this value. 10x is sufficient for most data; lower if error rate is very low."
    }

    # Estimate disk size required
    Int input_read_size  = ceil(size(input_reads, "GB"))
    Int final_disk_dize   = input_read_size * 5 + addldisk

    # Create output file name
    String output_counts_fn = "${sample_id}_${read_type}_counts.txt"

    command <<<
        set -eux -o pipefail

        READS=(~{sep=" " input_reads})

        ## dump converted reads into output folder
        mkdir output

        ## convert all reads to fastq.gz format
        for READFILE in "${READS[@]}"
            do
                echo "$READFILE"

                FILENAME="$(basename -- $READFILE)"
                PREFIX="${FILENAME%.*}"
                SUFFIX="${FILENAME##*.}"

                if [[ "$SUFFIX" == "bam" ]] ; then
                    samtools fastq -@~{threadCount} $READFILE | gzip --fast --stdout > output/${PREFIX}.fastq.gz

                elif [[ "$SUFFIX" == "cram" ]] ; then
                    if [[ ! -f "~{cram_reference}" ]] ; then
                        echo "Could not extract $FILENAME, reference file not supplied"
                        exit 1
                    fi
                    ln -s ~{cram_reference}
                    samtools fastq -@~{threadCount} --reference `basename ~{cram_reference}` $READFILE | gzip --fast --stdout >  output/${PREFIX}.fastq.gz
                elif [[ "$SUFFIX" == "fastq" ]] ; then
                    gzip $READFILE > output/${PREFIX}.fastq.gz

                elif [[ "$SUFFIX" != "gz" ]] ; then
                    echo "Unsupported file type: ${PREFIX}.${SUFFIX}"
                    exit 1
                else
                    ln -s $READFILE output/$FILENAME
                fi

            done


        ## Count number of k-mers found
        /opt/ntsm/bin/ntsmCount \
            -s /opt/ntsm/data/human_sites_n10.fa \
            -t ~{threadCount} \
            -m ~{max_coverage} \
            output/*.fastq.gz \
            > ~{output_counts_fn}
    >>>

    output {
        File ntsm_counts = output_counts_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + final_disk_dize + " SSD"
        docker: "iviolich/ntsm:1.2.1"
        preemptible: preempts
    }
}

task ntsm_eval {
    input {
        Array[File] count_files
        String output_prefix

        Int memSizeGB   = 4
        Int threadCount = 4
        Int diskSize    = 32
        Int preempts    = 2
    }

    parameter_meta {
        count_files: "Count files produced by ntsm_count, one per sample/read-type. All files are evaluated together."
        output_prefix: "Prefix for the output eval file."
    }

    String output_eval_fn = "${output_prefix}_ntsm_eval.tsv"

    command <<<
        set -eux -o pipefail

        ## Symlink count files by basename so ntsmEval uses clean sample names
        COUNT_FILES=(~{sep=" " count_files})
        for f in "${COUNT_FILES[@]}"; do
            ln -s "$f" "$(basename $f)"
        done

        /opt/ntsm/bin/ntsmEval \
            --all \
            $(for f in "${COUNT_FILES[@]}"; do basename "$f"; done) \
            > ~{output_eval_fn}
    >>>

    output {
        File ntsm_eval_out = output_eval_fn
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        docker: "iviolich/ntsm:1.2.1"
        preemptible: preempts
    }
}
