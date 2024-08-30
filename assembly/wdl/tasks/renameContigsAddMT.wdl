version 1.0

workflow renameContigsAddMT_wf {
    
    meta {
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Workflow to add mito contig and rename fasta headers for Genbank."
    }

    parameter_meta {
        sampleName: "Sample or assembly name. Used for naming output."
        outputFileTag: "Identifier to add in addition to sample name. Used for naming output."
        haplotype: "Haplotype of assembly. Allowed values are 1/2. Mito is only added to haplotype 2."
        inputFastaGZ: "Assembly to rename and add mito to. Should have been previously stripped of all mito. Input must be fasta that is gzipped."
        t2t_sequences: "Sequence IDs which represent an chromsome T2T. Must have header and column 1 is sequence IDs, column 4 is chromosome assignment formatted as chr1, for example"
        mitoAssembly: "Unzipped fasta assembly of mitochondrion."
    }
    
    input {
        String sampleName
        String outputFileTag
        Int haplotype
        File inputFastaGZ
        File t2t_sequences

        File? mitoAssembly
    }

    call renameContigsAddMT {
        input:
            sampleName = sampleName,
            outputFileTag = outputFileTag,
            haplotype = haplotype,
            inputFastaGZ = inputFastaGZ,
            t2t_sequences = t2t_sequences,
            mitoAssembly = mitoAssembly
    }

    output {
        File FinalAssembly  = renameContigsAddMT.FinalAssembly
    }
}

task renameContigsAddMT {

    input {
        String sampleName
        String outputFileTag
        Int haplotype
        File inputFastaGZ
        File t2t_sequences

        File? mitoAssembly

        Int threadCount    = 2
        Int memSizeGB      = 4
        Int diskSizeGB     = 64
        String dockerImage = "humanpangenomics/mitohifi_utils@sha256:cbd31256f13772d442f80a2e37cd2f36703a780a852d59dd6028e2cf8892c7d7"
    }


    String unzippedOrigFa  = basename(inputFastaGZ, ".gz")
    String outputFastaGz   = "~{sampleName}.~{outputFileTag}.fa.gz"
    String outputDupIDFile = "~{sampleName}_~{outputFileTag}_dup_ids.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{unzippedOrigFa}

        ## Check if we have a mito assembly (not all samples do). If so, add it to maternal haplotype
        if [ -s ~{mitoAssembly} ]; then

            ## Only add to haplotype 2
            ## By convention hap1 is paternal in trio, has chrY for HiC phased male samples
            ## hap2 is maternal in trio, has chrM, has chrX for HiC phased male samples
            if [[ ~{haplotype} == 2 ]]; then

                ## Rename just the mito contig (Genbank requires that it is labeled)
                # For example:   ">ptg000005l_rotated ptg000005l"
                ## Sould become: ">ptg000005l_rotated [location=mitochondrion]"
                sed 's/^\(>[^ ]*\).*/\1 [location=mitochondrion]/' ~{mitoAssembly} > renamed_mito.fasta 

                ## Add renamed mito contig to original assembly
                cat ~{unzippedOrigFa} renamed_mito.fasta > mito_added.fasta
            else
                ln -s ~{unzippedOrigFa} mito_added.fasta
            fi
        else
            ## There is no mito assembly, just link the file
            ln -s ~{unzippedOrigFa} mito_added.fasta
        fi 


        ## annotate fasta header so sequences that represent an entire chromosome are
        ## flagged by Genbank...

        ## Input file is of the form:
        # name    length  Ns  chromosome
        # h2tg000005l 253066670   27673   chr1
        # h2tg000014l 242466534   0   chr2
        # h2tg000008l 200287278   0   chr3

        ## Convert into:
        # h2tg000005l [chromosome=1]
        # h2tg000014l [chromosome=2]
        awk 'BEGIN {OFS="\t"} NR > 1 {gsub("chr", "", $4); print $1, "[location=chromosome][chromosome=" $4 "]"}' \
            ~{t2t_sequences} \
            > contig_to_chrom_map.txt

        ## Now actually replace the fasta headers...
        seqkit replace -p '(.+)' -r '$1 {kv}' -k mapping.txt \
            mito_added.fasta \
            --kv-file contig_to_chrom_map.txt \
            > mito_added_renamed_header.fasta
        

        ## Rename contig names to sampleName#1/2#contigName format (1 = paternal, 2 = maternal)
        sed "s/^>/>~{sampleName}\#~{haplotype}\#/" mito_added_renamed_header.fasta \
            > mito_added_renamed_header_sample_named.fasta

        cat mito_added_renamed_header_sample_named \
            | seqkit rmdup \
                --by-seq \
                -o ~{outputFastaGz} \
                --dup-num-file ~{outputDupIDFile}
    >>>

    output {

        File FinalAssembly  = outputFastaGz
        File dupIDFile      = outputDupIDFile
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}