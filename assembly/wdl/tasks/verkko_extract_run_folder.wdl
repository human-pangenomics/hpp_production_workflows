version 1.0


workflow verkko_extract_run_folder {

    call extract_run_folder
    
    output{
        File haplotype1_fa = extract_run_folder.haplotype1_fa
        File haplotype2_fa = extract_run_folder.haplotype2_fa
        File unassigned_fa = extract_run_folder.unassigned_fa
    }
}

task extract_run_folder {
    
    input {
        File run_folder
        String sample_name = "assembly"
        String tag         = "verkko_asm"

        Int threadCount = 6
        Int memSizeGB   = 24
        Int diskSizeGB  = 1000
        Int preemptible = 1
    }

    command <<<

        set -eux -o pipefail

        ## extract tar w/ snakemake run to cwd
        tar xvf ~{run_folder} --directory ./

        ## Get the name of the assembly folder
        asm_folder=$(ls -d */)

        ## extract haplotypes
        # awk < ${asm_folder}/final_*/assembly.fasta \
        #    'BEGIN {{
        #       FS="[ \\t]+"; OFS="\\t";
        #     }}
        #     (\$1 ~ /^haplotype1[0-9]+$/)    {{ print \$2 > "./haplotype1.names"; }}
        #     (\$2 ~ /^haplotype2[0-9]+$/)    {{ print \$2 > "./haplotype2.names"; }}
        #     (\$2 ~ /^unassigned-[0-9]+$/)   {{ print \$2 > "./unassigned.names"; }}'

        grep "haplotype1.*" ${asm_folder}/final_*/assembly.fasta > haplotype1.names
        grep "haplotype2.*" ${asm_folder}/final_*/assembly.fasta > haplotype2.names
        grep "unassigned.*" ${asm_folder}/final_*/assembly.fasta > unassigned.names

        seqtk subseq ${final_asm_folder}/assembly.fasta haplotype1.names | gzip > ~{sample_name}_~{tag}_haplotype1.fasta.gz &
        seqtk subseq ${final_asm_folder}/assembly.fasta haplotype2.names | gzip > ~{sample_name}_~{tag}_haplotype2.fasta.gz &
        seqtk subseq ${final_asm_folder}/assembly.fasta unassigned.names | gzip > ~{sample_name}_~{tag}_unassigned.fasta.gz &

        wait

    >>>

    output {

        File haplotype1_fa = "~{sample_name}_~{tag}_haplotype1.fasta.gz"
        File haplotype2_fa = "~{sample_name}_~{tag}_haplotype2.fasta.gz"
        File unassigned_fa = "~{sample_name}_~{tag}_unassigned.fasta.gz"
    }

    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: "quay.io/biocontainers/seqtk:1.3--hed695b0_2"
        preemptible: preemptible
    }
}