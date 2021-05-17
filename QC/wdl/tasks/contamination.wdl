version 1.0


#workflow runContamination {
#
#    input {
#        File assemblyFasta
#        File vecscreenContaminationDatabase
#        String dockerImage="tpesout/hpp_blast:latest"
#    }
#
#    call contaminationVecscreen {
#        input:
#            assemblyFasta=assemblyFasta,
#            vecscreenContaminationDatabase=vecscreenContaminationDatabase,
#            dockerImage=dockerImage
#    }
#
#	output {
#	    File vecscreenContamination = contaminationVecscreen.outputVecscreen
#	}
#}


workflow runContamination {

    input {
        String assemblyIdentifier
        File assemblyFasta
        File eukContaminationDatabase
        File mitoContaminationDatabase
        File plastidsContaminationDatabase
        File vecscreenContaminationDatabase
        File rrnaContaminationDatabase
        Array[File] refseqContaminationDatabases
        String dockerImage="tpesout/hpp_blast:latest"
    }

    call contaminationEuk {
        input:
            assemblyFasta=assemblyFasta,
            eukContaminationDatabase=eukContaminationDatabase,
            dockerImage=dockerImage
    }
    call contaminationMito {
        input:
            assemblyFasta=assemblyFasta,
            mitoContaminationDatabase=mitoContaminationDatabase,
            dockerImage=dockerImage
    }
    call contaminationPlastids {
        input:
            assemblyFasta=assemblyFasta,
            plastidsContaminationDatabase=plastidsContaminationDatabase,
            dockerImage=dockerImage
    }
    call contaminationVecscreen {
        input:
            assemblyFasta=assemblyFasta,
            vecscreenContaminationDatabase=vecscreenContaminationDatabase,
            dockerImage=dockerImage
    }
    call contaminationRRNA {
        input:
            assemblyFasta=assemblyFasta,
            rrnaContaminationDatabase=rrnaContaminationDatabase,
            dockerImage=dockerImage
    }
    call contaminationWindowmasker {
        input:
            assemblyFasta=assemblyFasta,
            dockerImage=dockerImage
    }
    scatter (refseqContaminationDatabase in refseqContaminationDatabases) {
        call contaminationRefseq {
            input:
                assemblyWindowmaskedFasta=contaminationWindowmasker.outputWindowmasker,
                refseqContaminationDatabase=refseqContaminationDatabase,
                dockerImage=dockerImage
        }
    }
    call mergeContaminationResults {
        input:
            assemblyIdentifier=assemblyIdentifier,
            eukOut=contaminationEuk.outputEuk,
            mitoOut=contaminationMito.outputMito,
            plastidsOut=contaminationPlastids.outputPlastids,
            rrnaOut=contaminationRRNA.outputRRNA,
            refseqOuts=contaminationRefseq.outputRefseq,
            vecscreenOut=contaminationVecscreen.outputVecscreen,
            dockerImage=dockerImage
    }

	output {
	    File eukContamination = contaminationEuk.outputEuk
	    File mitoContamination = contaminationMito.outputMito
	    File plastidsContamination = contaminationPlastids.outputPlastids
	    File windowmaskedFasta=contaminationWindowmasker.outputWindowmasker
	    Array[File] refseqContamination = contaminationRefseq.outputRefseq
	    File mergedResult = mergeContaminationResults.outputSummary
	}
}


task contaminationEuk {
    input {
        File assemblyFasta
        File eukContaminationDatabase
        String eukExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{eukContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{eukContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{eukContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -outfmt '7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore' \
            -perc_identity 90.0  ~{eukExtraArguments} \
            | awk '($3>=98.0 && $4>=50)||($3>=94.0 && $4>=100)||($3>=90.0 && $4>=200)' \
            > $PREFIX.contam_in_euks.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore\tsubject description" >$PREFIX.contam_in_euks.tsv
        cat $PREFIX.contam_in_euks.awk | { grep -v "^#" || true; } >>$PREFIX.contam_in_euks.tsv

	>>>
	output {
		File outputEuk = glob("*.contam_in_euks.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationMito {
    input {
        File assemblyFasta
        File mitoContaminationDatabase
        String mitoExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{mitoContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{mitoContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{mitoContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -perc_identity 90 \
            -soft_masking true \
            -outfmt 7 ~{mitoExtraArguments} \
            | awk '$4>=500' \
            > $PREFIX.mito.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore\tsubject description" >$PREFIX.mito.tsv
        cat $PREFIX.mito.awk | { grep -v "^#" || true; } >>$PREFIX.mito.tsv

	>>>
	output {
		File outputMito = glob("*.mito.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationPlastids {
    input {
        File assemblyFasta
        File plastidsContaminationDatabase
        String plastidsExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{plastidsContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{plastidsContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{plastidsContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -perc_identity 90 \
            -soft_masking true \
            -outfmt 7 ~{plastidsExtraArguments} \
            | awk '$4>=500' \
            > $PREFIX.plastids.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore\tsubject description" >$PREFIX.plastids.tsv
        cat $PREFIX.plastids.awk | { grep -v "^#" || true; } >>$PREFIX.plastids.tsv

	>>>
	output {
		File outputPlastids = glob("*.plastids.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationRRNA {
    input {
        File assemblyFasta
        File rrnaContaminationDatabase
        String rrnaExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{rrnaContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{rrnaContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{rrnaContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -template_length 18 \
            -template_type coding \
            -window_size 120 \
            -word_size 12 \
            -xdrop_gap 20 \
            -no_greedy \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 1E-9 \
            -gapextend 2 \
            -gapopen 4 \
            -penalty -5 \
            -perc_identity 93 \
            -reward 2 \
            -soft_masking true \
            -outfmt 7 \
            ~{rrnaExtraArguments} \
            | awk '$4>=100' \
            > $PREFIX.rrna.awk

        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore\tsubject description" >$PREFIX.rrna.tsv
        cat $PREFIX.rrna.awk | { grep -v "^#" || true; } >>$PREFIX.rrna.tsv

	>>>
	output {
		File outputRRNA = glob("*.rrna.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task contaminationWindowmasker {
    input {
        File assemblyFasta
        String refseqExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        windowmasker -mk_counts -in $ASM_FILENAME -out $PREFIX.stage1
        windowmasker -ustat $PREFIX.stage1 -in $ASM_FILENAME -out $PREFIX.masked.fa -dust true -outfmt fasta

	>>>
	output {
		File outputWindowmasker = glob("*.masked.fa")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}

task contaminationRefseq {
    input {
        File assemblyWindowmaskedFasta
        File refseqContaminationDatabase
        String refseqExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 128
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyWindowmaskedFasta}")
        ln -s ~{assemblyWindowmaskedFasta}
        PREFIX="${ASM_FILENAME%.masked.fa}"
        RSDB_ID=$(basename -- "~{refseqContaminationDatabase}")
        RSDB_ID="${RSDB_ID%.*}"

        # make db index
        DB_FILENAME=$(basename -- "~{refseqContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{refseqContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{refseqContaminationDatabase}
        fi
        makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # blastoff
        blastn -query $ASM_FILENAME \
            -db $DB_FILENAME \
            -task megablast \
            -num_threads ~{threadCount} \
            -word_size 28 \
            -best_hit_overhang 0.1 \
            -best_hit_score_edge 0.1 \
            -dust yes \
            -evalue 0.0001 \
            -min_raw_gapped_score 100 \
            -penalty -5 \
            -perc_identity 98 \
            -soft_masking true \
            -outfmt "7 std staxids sscinames" ~{refseqExtraArguments} \
            > $PREFIX.refseq.$RSDB_ID.out
#            -negative_taxidlist ${negative_tax_id_list} \

        echo -e "========== Vecscreen: $RSDB_ID ==========" >$PREFIX.refseq.$RSDB_ID.tsv
        echo -e "#query acc.ver\tsubject acc.ver\t% id\tlength\tmsmatch\tgaps\tq start\tq end\ts start\ts end\tevalue\tscore\tsubject description" >>$PREFIX.refseq.$RSDB_ID.tsv
        cat $PREFIX.refseq.$RSDB_ID.out | { grep -v "^#" || true; } >>$PREFIX.refseq.$RSDB_ID.tsv

	>>>
	output {
		File outputRefseq = glob("*.refseq.*.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}



task contaminationVecscreen {
    input {
        File assemblyFasta
        File vecscreenContaminationDatabase
        String vecscreenExtraArguments=""
        Int memSizeGB = 8
        Int threadCount = 8
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # chunk it (vecscreen cannot handle long sequences)
        chunk_assembly.py $ASM_FILENAME $ASM_FILENAME.chunked

        # index it (old version of blast.. don't know if we NEED 2.7.1 but the latest doesn't work)
        DB_FILENAME=$(basename -- "~{vecscreenContaminationDatabase}")
        if [[ $DB_FILENAME =~ \.gz$ ]]; then
            cp ~{vecscreenContaminationDatabase} .
            gunzip $DB_FILENAME
            DB_FILENAME="${DB_FILENAME%.gz}"
        else
            ln -s ~{vecscreenContaminationDatabase}
        fi
        /opt/blast/ncbi-blast-2.7.1+/bin/makeblastdb -in $DB_FILENAME -input_type fasta -dbtype nucl

        # screen it
        vecscreen -d $DB_FILENAME \
            -f3 \
            -i $ASM_FILENAME.chunked \
            -o $PREFIX.vecscreen.out ~{vecscreenExtraArguments}
        VSlistTo1HitPerLine.awk \
            suspect=0 \
            weak=0 \
            $PREFIX.vecscreen.out \
            > $PREFIX.vecscreen.filtered.out

        # unchunk it
        unchunk_vecscreen.py $PREFIX.vecscreen.filtered.out $PREFIX.vecscreen.tsv

	>>>
	output {
		File outputVecscreen = glob("*.vecscreen.tsv")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}


task mergeContaminationResults {
    input {
        String assemblyIdentifier
        File eukOut
        File mitoOut
        File rrnaOut
        File plastidsOut
        Array[File] refseqOuts
        File rrnaOut
        File vecscreenOut
        Int memSizeGB = 2
        Int threadCount = 1
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_blast:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        OUT="~{assemblyIdentifier}.contamination.txt"

        echo "========== Eukaryote ==========" >>$OUT
        cat ~{eukOut} >>$OUT

        echo "" >>$OUT
        echo "========== Mitocondria ==========" >>$OUT
        cat ~{mitoOut} >>$OUT

        echo "" >>$OUT
        echo "========== Plastids ==========" >>$OUT
        cat ~{plastidsOut} >>$OUT

        echo "" >>$OUT
        echo "========== RRNA ==========" >>$OUT
        cat ~{rrnaOut} >>$OUT

        echo "" >>$OUT
        echo "========== Vecscreen ==========" >>$OUT
        cat ~{vecscreenOut} >>$OUT

        echo "" >>$OUT
        echo "========== RefSeq ==========" >>$OUT
        cat ~{sep=" " refseqOuts} >>$OUT

	>>>
	output {
		File outputSummary = glob("*.contamination.txt")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
