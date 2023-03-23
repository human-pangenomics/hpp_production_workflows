version 1.0

workflow runAssemblyCleanUp{
    call assemblyCleanUp
    output {
        File miscFilesTarGz = assemblyCleanUp.miscFilesTarGz
        File mitoFastaGz = assemblyCleanUp.mitoFastaGz
        File ebvFastaGz = assemblyCleanUp.ebvFastaGz
        File rdnaFastaGz = assemblyCleanUp.rdnaFastaGz
        File cleanedFastaGz = assemblyCleanUp.cleanedFastaGz
    }
}


task assemblyCleanUp {
    input{
        File paternalGfaTarGz
        File maternalGfaTarGz
        String childID
        # runtime configurations
        Int memSizeGB=32
        Int threadCount=8
        Int diskSizeGB=128
        Int preemptible=1
        String dockerImage="quay.io/biocontainers/verkko:1.3.1--h64afbab_0"
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

        mkdir outputs
        tar -xvzf ~{paternalGfaTarGz} --strip-components 1 --directory outputs
        tar -xvzf ~{maternalGfaTarGz} --strip-components 1 --directory outputs
        cd outputs 

        cat *.dip.hap[12].p_ctg.gfa | awk '{if (match($1, "^S")) { print ">"$2; print $3}}' > combined.fasta
        cat *.dip.hap[12].p_ctg.gfa > combined.gfa
        cat *.dip.hap[12].p_ctg.noseq.gfa | awk '{if (match($1, "^S")) { print "path "$2"\t"$2; print "piece000002"; print "end"; }}' > combined.scfmap
        cat *.dip.hap[12].p_ctg.noseq.gfa | awk -v NAME="" '{if (match($1, "^S")) { if (NAME != "") print NAME"\t"SUM/LEN; SUM=0; NAME=$2; LEN=substr($4, 6, length($4)); } if (match($1, "^A")) { SUM+=$7-$6}} END {print NAME"\t"SUM/LEN}' > combined.csv
      
         /usr/local/lib/verkko/scripts/screen-assembly.pl \
          --assembly      combined.fasta \
          --graph         combined.gfa \
          --graphmap      combined.scfmap \
          --hifi-coverage combined.csv \
          --minlength     100000 \
          --contaminant   ebv /usr/local/lib/verkko/data/human-ebv-AJ507799.2.fasta.gz mito /usr/local/lib/verkko/data/human-mito-NC_012920.1.fasta.gz rdna /usr/local/lib/verkko/data/human-rdna-KY962518.1.fasta.gz \
          --output        assembly \
          --mashmap       mashmap \
          --threads       ~{threadCount} \
          > ./screen-assembly.out 2> ./screen-assembly.err

          #  Try to circularize, the script will output to the same file as input since it caches the sequences
          #  if it fails to circularize it will output sequence unchanged
          for xx in assembly.*.exemplar.fasta ; do
              if [ -e $xx ]; then
                  python3 /usr/local/lib/verkko/scripts/circularize_ctgs.py -p 10 -f 0.90 -o $xx --min-ovl 2500 $xx
              fi
          done
          cd ../
         
          mv outputs/assembly.fasta ~{childID}.cleaned.fa
          gzip ~{childID}.cleaned.fa
 
          for suffix in mito rdna ebv; do
              mv outputs/assembly.${suffix}.fasta ~{childID}.${suffix}.fa
              gzip ~{childID}.${suffix}.fa 
          done

          mv outputs ~{childID}.clean_up
          tar -cf ~{childID}.clean_up.tar ~{childID}.clean_up
          gzip ~{childID}.clean_up.tar
          
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: preemptible
    }

    output {
        File miscFilesTarGz = "${childID}.clean_up.tar.gz"
        File mitoFastaGz = "${childID}.mito.fa.gz"
        File ebvFastaGz = "${childID}.ebv.fa.gz"
        File rdnaFastaGz = "${childID}.rdna.fa.gz"
        File cleanedFastaGz = "${childID}.cleaned.fa.gz"
    }
}

