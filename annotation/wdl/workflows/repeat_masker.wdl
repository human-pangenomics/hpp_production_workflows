version 1.0
# WDL ReapeatMasker workflow - runs each contig/chromosome individually 

import "../tasks/bed_to_bigbed.wdl" as bigbed_wf

workflow RepeatMasker {
    input {
         File fasta
         
         File? additionalRMModels
         File? rmskAutoSql
         File? rmskAlignAutoSql
         String fName = sub(basename(fasta), "\.(fa|fasta)(\.gz)?$", "")
     }
    
    call createArray {
        input:
            fasta=fasta
    }

    scatter (subFasta in createArray.contigArray) {
        call maskContig {
            input:
                subsequence=subFasta,
                additionalRMModels=additionalRMModels
        }

        call outToBed {
            input:
                maskedOut = maskContig.outFile
        }
    }

    call finalizeFiles {
        input:
            bedFiles       = outToBed.RMbed,
            maskedFastas   = maskContig.maskedFa,
            subseqTarGZs   = maskContig.subseqTarGZ,
            fName          = fName
    }

    
    if (defined(rmskAutoSql) && defined(rmskAlignAutoSql)) {

        call bigbed_wf.bed_to_bigbed_wf as bed_to_bigbed_rmsk {
            input:
                input_bed   = finalizeFiles.rmskBed,
                auto_sql    = select_first([rmskAutoSql]),
                assembly_fa = fasta,
                type_str    = "bed9+5",
                tab_delim   = true 
        }

        call bigbed_wf.bed_to_bigbed_wf as bed_to_bigbed_align_rmsk {
            input:
                input_bed   = finalizeFiles.rmskAlignBed,
                auto_sql    = select_first([rmskAlignAutoSql]),
                assembly_fa = fasta,
                type_str    = "bed3+14",
                tab_delim   = true 
        }        
    }

    output {
        File repeatMaskerBed        = finalizeFiles.repeatMaskerBed
        File rmskBed                = finalizeFiles.rmskBed
        File rmskAlignBed           = finalizeFiles.rmskAlignBed
        File finalMaskedFasta       = finalizeFiles.finalMaskedFasta        
        File repeatMaskerOutFile    = finalizeFiles.repeatMaskerOutFile
        File repeatMaskerTarGZ      = finalizeFiles.repeatMaskerTarGZ

        File? rmskBigBed            = bed_to_bigbed_rmsk.output_bigbed
        File? rmskAlignBigBed       = bed_to_bigbed_align_rmsk.output_bigbed
    }

    parameter_meta {
         fasta: "Assembly for annotation. (supported formats: .fasta, .fa, .fasta.gz, .fa.gz)"
         additionalRMModels: "(Optional) Models to add to Dfam DB before masking."
         rmskAutoSql: "(Optional) AutoSql file for converting rmsk file to bigbed. See https://genome.ucsc.edu/goldenPath/help/bigRmsk.html"
         rmskAlignAutoSql: "(Optional) AutoSql file for converting rmskAlign file to bigbed. See https://genome.ucsc.edu/goldenPath/help/bigRmsk.html"
    }

    meta {
        author: "Hailey Loucks"
        email: "hloucks@ucsc.edu"
    }
}

task createArray {
    input {
        File fasta
	String fName

        Int threadCount = 2
        Int memSizeGB   = 16
        Int diskSize    = 32
        Int preemptible = 1
    }
    command <<<

        ## first check if assembly_fa needs to be unzipped 
        if [[ ~{fasta} =~ \.gz$ ]]; then
            gunzip -fc $fasta > ~{fName}.fa
        else
            cat ~{fasta} > ~{fName}.fa
        fi 

        ## split fasta and name outputs with sequence ID.
        ## Note that "#" are converted to "_" in the file names (but not the headers) 
        ## to prevent the "#" being interpreted URL encoded (%23)
        awk '/^>/ { file=substr($1,2); gsub("#","_",file); file=file ".fa" } { print > file }' ~{fName}.fa

    >>>
    output {
        Array[File] contigArray = glob("*.fa")
    }
    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSize + " SSD"
        docker: "ubuntu@sha256:152dc042452c496007f07ca9127571cb9c29697f42acbfad72324b2bb2e43c98" # 18.04
        preemptible : preemptible
    }
}


task maskContig {
    input {
        File subsequence
        File? additionalRMModels
        String subsequenceName=basename(subsequence)

        # 4 threads per -pa 1 in RepeatMasker call required, recommend at least 32 threads with RM command at -pa 8
        Int threadCount = 8
        Int memSizeGB   = 16
        Int diskSize    = 32
        Int preemptible = 1
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        # each parallel job uses 4 threads, budget 4 threads per job, rounding up
        NPARALLEL=$(( (~{threadCount} + 3) / 4 ))
        
        ## If additional models are passed, add them to the Dfam DB
        ## For example: Add new models from Primate X/Y and Human Y papers 
        ## to Dfam 3.7 (which is included in dfam/tetools:1.85)
        if [ -n "~{additionalRMModels}" ]; then
            wd=$(pwd)
            cd /opt/RepeatMasker/Libraries/
            python3 ../famdb.py -i ./RepeatMaskerLib.h5 append ~{additionalRMModels} --name 'homo_sapiens'
            cd $wd
        fi 

        ## -xsmall: output soft mask
        ## -align: output align files (for big rmsk files)
        RepeatMasker -s \
            -pa ${NPARALLEL} \
            -e ncbi \
            ~{subsequence} \
            -species human \
            -xsmall \
            -libdir /opt/RepeatMasker/Libraries/ \
            -align \
            -dir .

        # for small contigs - if there are no repeats found put the unmasked 
        ## sequence in and create empty files for rest of output
        if ! test -f ~{subsequenceName}.masked; then 
            touch ~{subsequenceName}.out \
            && touch ~{subsequenceName}.tbl \
            && touch ~{subsequenceName}.align \
            && cat ~{subsequence} > ~{subsequenceName}.masked \
            && touch ~{subsequenceName}.rmsk.bed \
            && touch ~{subsequenceName}.rmsk.align.bed
        else
            ## Create rmsk files for better genome browser viewing
            ## See https://genome.ucsc.edu/goldenPath/help/bigRmsk.html
            /opt/RepeatMasker/util/rmToTrackHub.pl \
                -out ~{subsequenceName}.out \
                -align ~{subsequenceName}.align 

            ## rename to reflect bed file format
            mv ~{subsequenceName}.join.tsv  ~{subsequenceName}.rmsk.bed
            mv ~{subsequenceName}.align.tsv ~{subsequenceName}.rmsk.align.bed
        fi

        ## create archive of all outputs for portability
        mkdir ~{subsequenceName}
        
        cp ~{subsequenceName}.{out,tbl,align,masked,rmsk.bed,rmsk.align.bed} ~{subsequenceName}/

        tar -czf ~{subsequenceName}.tar.gz ~{subsequenceName}
        
    >>>

    output {
         File outFile      = "~{subsequenceName}.out"
         File tblFile      = "~{subsequenceName}.tbl"
         File alignFile    = "~{subsequenceName}.align"
         File maskedFa     = "~{subsequenceName}.masked"
         File rmskBed      = "~{subsequenceName}.rmsk.bed"
         File rmskAlignBed = "~{subsequenceName}.rmsk.align.bed"
         File subseqTarGZ  = "~{subsequenceName}.tar.gz"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSize + " SSD"
        docker: 'humanpangenomics/dfam_tetools@sha256:31fa91da05360fbf1b8a7fa4f011093b3b1771a8a338e79b93fc1a6777d01781' # 1.85
        preemptible : preemptible
    }
}

task outToBed {
    input {
        File maskedOut
        String subsequenceName=basename(maskedOut, ".out")

        Int threadCount = 2
        Int memSizeGB   = 16
        Int diskSize    = 32
        Int preemptible = 1
    }
    command <<<

        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace
        
        /opt/RM2Bed.py \
            ~{maskedOut} \
            --out_dir . \
            --out_prefix ~{subsequenceName} \
            --ovlp_resolution 'higher_score'

    >>>
    output {
        File RMbed = "~{subsequenceName}_rm.bed"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        disks: "local-disk " + diskSize + " SSD"
        docker: 'humanpangenomics/rm2bed@sha256:a5d415c3aa6372c9d3f725a11538d0e1ef411fb0810cf0c063346a0e4be5b4a0'
        preemptible : preemptible
    }
}

task finalizeFiles {
    input {
        Array[File] bedFiles
        Array[File] maskedFastas
        Array[File] subseqTarGZs
        String fName

        Int threadCount = 2
        Int memSizeGB   = 16
        Int diskSize    = 96
        Int preemptible = 1
    }
    command <<<
        #handle potential errors and quit early
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # Extract all archives
        mkdir tar_directory
        
        for tarfile in ~{sep=' ' subseqTarGZs}; do
            tar xzf ${tarfile} -C tar_directory
        done

        #concatenate the .out files 
        find ./tar_directory -name "*.out" | xargs cat > rm.tmp
        head -n 3 rm.tmp > ~{fName}_repeat_masker.out
        find ./tar_directory -name "*.out" -exec sed '1,3d' {} \; >> ~{fName}_repeat_masker.out

        # concatenate the rmsk bed files
        find ./tar_directory -name "*.rmsk.bed" | xargs cat > ~{fName}.rmsk.bed
        bedtools sort -i ~{fName}.rmsk.bed > ~{fName}_repeat_masker_rmsk.bed

        # concatenate the rmsk align bed files
        find ./tar_directory -name "*.rmsk.align.bed" | xargs cat > ~{fName}.rmsk.align.bed
        bedtools sort -i ~{fName}.rmsk.align.bed > ~{fName}_repeat_masker_rmsk.align.bed


        # concatenate the masked fastas
        cat ~{sep=' ' maskedFastas} | seqkit sort --natural-order --two-pass | pigz > ~{fName}_repeat_masker_masked.fasta.gz


        # concatenate the RM2BED bed files
        cat ~{sep=' ' bedFiles} > ~{fName}.bed
        bedtools sort -i ~{fName}.bed > ~{fName}_repeat_masker.bed


        # Create an archive of the sequence-level tbl, out, and align files
        mkdir -p ~{fName}_repeat_masker/{out,tbl,align}
        
        find ./tar_directory -name "*.out"   -exec cp {} ~{fName}_repeat_masker/out/ \;
        find ./tar_directory -name "*.tbl"   -exec cp {} ~{fName}_repeat_masker/tbl/ \;
        find ./tar_directory -name "*.align" -exec cp {} ~{fName}_repeat_masker/align/ \;

        tar -czf ~{fName}_repeat_masker.tar.gz ~{fName}_repeat_masker

    >>> 
    output {
        File repeatMaskerBed     = "~{fName}_repeat_masker.bed"
        File repeatMaskerOutFile = "~{fName}_repeat_masker.out"        
        File rmskBed             = "~{fName}_repeat_masker_rmsk.bed"
        File rmskAlignBed        = "~{fName}_repeat_masker_rmsk.align.bed"
        File finalMaskedFasta    = "~{fName}_repeat_masker_masked.fasta.gz"
        File repeatMaskerTarGZ   = "~{fName}_repeat_masker.tar.gz"
    }

    runtime {
        cpu: threadCount        
        memory: memSizeGB + " GB"
        preemptible : preemptible
        disks: "local-disk " + diskSize + " SSD"
        docker: "humanpangenomics/sequence_toolbox@sha256:07fd158dce6b25b1731631cec161e554d0704cbf1227edb890b778c6b6fcb39b"
    }
}
