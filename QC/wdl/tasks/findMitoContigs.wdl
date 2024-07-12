version 1.0

workflow findMitoContigs {
    meta{
        author: "Julian Lucas"
        email: "juklucas@ucsc.edu"
        description: "Workflow to find mitochondrial sequences in an assembly using BLAST. NUMTs are ignored. Process is based on MitoHiFi and NCBI's internal BLAST search implemented for Genbank Screen."
    }

    input {
        String sample_id
        String haplotype
        File inputFastaGZ
        File mito_ref
        File ncbi_mito_blast_db
    }
    
    parameter_meta {
        sample_id: "Sample or assembly name for fasta that is being searched. Used for output naming"
        haplotype: "Haplotype or additional info for fasta that is being searched. Used for output naming"
        inputFastaGZ: "Assembly to search for Mitochondrial sequences."
        mito_ref: "Mito fasta from same, or closely related, species. Cannot be gzipped. Used to make BLAST DB"
        ncbi_mito_blast_db: "NCBI Genomic Mitochondrial Reference Sequences as a BLAST DB. [See metadata](https://ftp.ncbi.nlm.nih.gov/blast/db/mito-nucl-metadata.json)"
    }

    call blastFasta {
        input:
            mito_ref     = mito_ref,
            inputFastaGZ = inputFastaGZ,
            sample_id    = sample_id,
            haplotype    = haplotype
    }

    call ncbiBlastFasta {
        input:
            ncbi_mito_blast_db = ncbi_mito_blast_db,
            inputFastaGZ       = inputFastaGZ,
            sample_id          = sample_id,
            haplotype          = haplotype
    }

    call parseBlastOutput as parseBlast {
        input:
            sample_id   = "~{sample_id}_~{haplotype}",
            tag         = "mito",
            blastOutput = blastFasta.blastOutput
    }

    call parseBlastOutput as parseNcbiBlast {
        input:
            sample_id   = "~{sample_id}_~{haplotype}",
            tag         = "ncbi",
            blastOutput = ncbiBlastFasta.blastOutput
    }

    call pullContigIDs {
        input:
            mito_hits  = [parseBlast.parsedBlastOutput, parseNcbiBlast.parsedBlastOutput],
            sample_id  = sample_id,
            tag        = haplotype
    }

    output {
        ## deduplicated contig list found with blast and NCBI's blast search
        File mitoContigs       = pullContigIDs.mitoContigs

        ## output from blast search for mito
        File blastOutput       = blastFasta.blastOutput
        File parsedBlast       = parseBlast.parsedBlastOutput

        ## output from an emulation of NCBI's blast search
        File ncbiBlastOutput   = ncbiBlastFasta.blastOutput
        File ncbiParsedBlast   = parseNcbiBlast.parsedBlastOutput
    }
}


task blastFasta {

    input {
        File mito_ref
        String sample_id
        String haplotype
        File inputFastaGZ
        
        Int threadCount    = 2
        Int memSizeGB      = 4
        Int diskSizeGB     = 64
        String dockerImage = "ncbi/blast@sha256:77a24a340683c2f4883e2d5295bf63277743579239ada939370c19ca5622ef5f" # 2.15.0
    }

    String inputFasta = basename(inputFastaGZ, ".gz")
    String blastOutputName = "${sample_id}.${haplotype}.mito_blast_out.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## Create blastn database
        makeblastdb \
            -in ~{mito_ref} \
            -title chrM_Reference \
            -dbtype nucl \
            -out chrM_ref/chrM_ref


        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{inputFasta}
        

        ## Query fasta against MT database
        blastn \
            -query ~{inputFasta} \
            -db chrM_ref/chrM_ref \
            -num_threads 2 \
            -dust no \
            -soft_masking false \
            -out ~{blastOutputName} \
            -outfmt '6 std qlen slen'

    >>>

    output {
        File blastOutput = blastOutputName
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task ncbiBlastFasta {

    input {
        File ncbi_mito_blast_db
        String sample_id
        String haplotype
        File inputFastaGZ
        
        Int threadCount    = 2
        Int memSizeGB      = 8
        Int diskSizeGB     = 64
        String dockerImage = "ncbi/blast@sha256:77a24a340683c2f4883e2d5295bf63277743579239ada939370c19ca5622ef5f" # 2.15.0
    }

    String inputFasta = basename(inputFastaGZ, ".gz")
    String blastOutputName = "${sample_id}.${haplotype}.ncbi_mito_blast_out.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## gunzip input fasta 
        gunzip -c ~{inputFastaGZ} > ~{inputFasta}
        

        mkdir -p mito_db
        tar -xzvf ~{ncbi_mito_blast_db} -C mito_db

        blastn \
            -query ~{inputFasta} \
            -num_threads 2 \
            -db mito_db/mito \
            -out % \
            -task megablast \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -perc_identity 98.6 \
            -soft_masking true \
            -outfmt '6 std qlen slen' \
            > ~{blastOutputName}

    >>>

    output {
        File blastOutput = blastOutputName
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task parseBlastOutput {

    input {
        String sample_id
        String tag
        File blastOutput
        Int threadCount    = 1
        Int memSizeGB      = 4
        Int diskSizeGB     = 64
        String dockerImage = "juklucas/parse_mito_blast:latest"
    }

    String parsedBlastOutputName = "${sample_id}.${tag}.parsedMitoBlast.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        ## parse blast output (aggregate blast hits into contig-level summary; 
        ## only outputs MT-only contigs). Slightly modified version of parse_blast.py 
        ## script found in MitoHiFi (https://github.com/marcelauliano/MitoHiFi)
        parse_mito_blast_results.py ~{blastOutput} ~{parsedBlastOutputName}
        
    >>>

    output {
        File parsedBlastOutput = parsedBlastOutputName
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task pullContigIDs {

    input {
        Array[File] mito_hits
        String sample_id
        String tag

        Int threadCount    = 1
        Int memSizeGB      = 4
        Int diskSizeGB     = 64
        String dockerImage = "juklucas/parse_mito_blast:latest"
    }

    String outputName = "${sample_id}.${tag}.mito_ids.txt"

    command <<<

        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # Create a temporary file to hold all first columns
        temp_file=$(mktemp)

        ## Expect input to be of the form:
        ## qseqid  %q_in_match     leng_query      s_length        perc
        ## h2tg000082l     100.06037553583289      16563   16569   99.96378779648741

        # Loop over all provided files and extract the first column
        for file in ~{sep=' ' mito_hits}; do
            ## skip header line, print just contig ID
            awk 'NR > 1 {print $1}' "${file}" >> ${temp_file}
        done

        # Sort and remove duplicates, then output the result
        sort ${temp_file} | uniq > ~{outputName}
        
    >>>

    output {
        File mitoContigs = outputName
    }

    runtime {
        cpu: threadCount
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

