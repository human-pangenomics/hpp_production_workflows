version 1.0

import "../tasks/ntsm.wdl" as ntsm_check

workflow illumina_wgs_qc_wf {
    input {
        File   wgs_reads
        String sample_id
        File?  cram_reference
    }

    parameter_meta {
        wgs_reads:      "Illumina WGS reads. Can be CRAM, BAM, or fastq.gz."
        sample_id:      "Sample ID used for output file naming."
        cram_reference: "Reference FASTA required if wgs_reads is a CRAM file."
    }

    meta {
        author: "Ivo Violich"
        email:  "iviolich@ucsc.edu"
    }

    ## count total bases and compute coverage
    call count_wgs_bases {
        input:
            wgs_reads       = wgs_reads,
            cram_reference  = cram_reference,
            sample_id       = sample_id
    }

    ## ntsm k-mer count for cross-sample swap detection
    call ntsm_check.ntsm_count as ntsm_wf {
        input:
            input_reads    = [wgs_reads],
            cram_reference = cram_reference,
            sample_id      = sample_id,
            read_type      = "wgs"
    }

    ## combine into summary TSV
    call summarize_wgs_qc {
        input:
            coverage_file = count_wgs_bases.coverage_file,
            file_name     = sample_id
    }

    output {
        File wgs_qc_summary  = summarize_wgs_qc.summary_file

        ## ntsm count file (batch eval run separately across all samples)
        File wgs_ntsm_counts = ntsm_wf.ntsm_counts
    }
}


task count_wgs_bases {
    input {
        File   wgs_reads
        String sample_id
        File?  cram_reference

        Int memSizeGB   = 8
        Int threadCount = 4
        Int addldisk    = 10
        Int preempts    = 2
    }

    Int input_size  = ceil(size(wgs_reads, "GB"))
    Int disk_size   = input_size * 3 + addldisk

    command <<<
        set -euo pipefail

        READS="~{wgs_reads}"
        SUFFIX="${READS##*.}"

        if [[ "$SUFFIX" == "cram" || "$SUFFIX" == "bam" ]]; then
            if [[ "$SUFFIX" == "cram" ]]; then
                if [[ ! -f "~{default="" cram_reference}" ]]; then
                    echo "ERROR: CRAM reference not supplied"
                    exit 1
                fi
                ln -s ~{cram_reference} .
                REF_ARG="--reference $(basename ~{cram_reference})"
            else
                REF_ARG=""
            fi
            TOTAL_BASES=$(samtools fastq -@ ~{threadCount} $REF_ARG "$READS" \
                | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
        else
            ## fastq or fastq.gz
            if [[ "$SUFFIX" == "gz" ]]; then
                TOTAL_BASES=$(zcat "$READS" | awk 'NR%4==2 {sum+=length($0)} END {print sum}')
            else
                TOTAL_BASES=$(awk 'NR%4==2 {sum+=length($0)} END {print sum}' "$READS")
            fi
        fi

        echo "$TOTAL_BASES" > ~{sample_id}_base_count.txt
    >>>

    output {
        File coverage_file = "~{sample_id}_base_count.txt"
    }

    runtime {
        memory:      memSizeGB + " GB"
        cpu:         threadCount
        disks:       "local-disk " + disk_size + " SSD"
        docker:      "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
        preemptible: preempts
    }
}


task summarize_wgs_qc {
    input {
        File   coverage_file
        String file_name

        Int memSizeGB   = 4
        Int threadCount = 1
        Int disk_size   = 16
        Int preempts    = 2
    }

    command <<<
        set -euo pipefail

        echo -e "file_name\ttotal_bp\tcoverage" > "~{file_name}.summary.tsv"

        cov_line=$(awk '
            BEGIN {
                genome_size = 3.3e9;
                OFS = "\t"
            }
            {
                total_bases = $1;
                coverage = total_bases / genome_size;
                printf("%d\t%.2f\n", total_bases, coverage);
            }
        ' "~{coverage_file}")

        echo -e "~{file_name}\t$cov_line" >> "~{file_name}.summary.tsv"
    >>>

    output {
        File summary_file = "~{file_name}.summary.tsv"
    }

    runtime {
        memory:      memSizeGB + " GB"
        cpu:         threadCount
        disks:       "local-disk " + disk_size + " SSD"
        docker:      "registry.access.redhat.com/ubi7/ubi:latest"
        preemptible: preempts
    }
}
