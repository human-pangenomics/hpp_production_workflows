version 1.0 

workflow runMapBlocks{
    call mapBlocks
    output {
        File mergedMappedBlocksBed = mapBlocks.mergedMappedBlocksBed
        File mappedBlocksBed = mapBlocks.mappedBlocksBed
        #File mappedLowMQBlocksBed = mapBlocks.mappedLowMQBlocksBed
        File unmappedBlocksBed = mapBlocks.unmappedBlocksBed
        #File skippedBlocksBed = mapBlocks.skippedBlocksBed
        File mergedRefBlocksBed = mapBlocks.mergedRefBlocksBed
        File refBlocksBed = mapBlocks.refBlocksBed
    }
}
task mapBlocks {
    input {
        File blocksBed
        File alignmentBam
        String suffix
        # runtime configurations
        Int memSize=8
        Int threadCount=2
        Int diskSize=32
        String dockerImage="quay.io/masri2019/hpp_asset:latest"
        Int preemptible=2
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
        
        FILENAME=$(basename ~{blocksBed})
        PREFIX=${FILENAME%.bed}
        cat ~{blocksBed} | awk '($3-$2)>10' > ${PREFIX}_gt10.bed
        samtools view -F256 -F4 -q20 ~{alignmentBam} | cut -f1-6 > no_seq.sam
        #samtools view -F256 -F4 ~{alignmentBam} | awk '$5 < 20' |  cut -f1-6 > mq_less_20_no_seq.sam
        samtools view -f4 ~{alignmentBam} | cut -f1 > unmapped_contigs.txt
        mkdir output
        python3 $MAP_BLOCKS_PY --sam no_seq.sam --bed ${PREFIX}_gt10.bed --outputContig output/contig.bed --outputMapped output/ref.bed --outputSkipped output/skipped.bed
        #python3 $MAP_BLOCKS_PY --sam mq_less_20_no_seq.sam --bed ~{blocksBed} --outputContig output/contig_20.bed --outputMapped output/ref_20.bed --outputSkipped output/skipped_20.bed
        #cat output/contig.bed output/contig_20.bed | bedtools sort -i - | bedtools merge -i - > output/contig_all.bed
        #cat output/skipped.bed output/skipped_20.bed | bedtools sort -i - | bedtools merge -i - > output/skipped_all.bed
        #bedtools subtract -a output/skipped_all.bed -b output/contig_all.bed | bedtools sort -i - | bedtools merge -i - > output/${PREFIX}.~{suffix}.skipped.bed
        #bedtools subtract -a output/contig_20.bed -b output/contig.bed | bedtools sort -i - | bedtools merge -i - > output/${PREFIX}.~{suffix}.mq_lq20.bed
        cat ${PREFIX}_gt10.bed | grep -f unmapped_contigs.txt | bedtools sort -i - | bedtools merge -i - > output/${PREFIX}.~{suffix}.unmapped.bed
        bedtools sort -i output/contig.bed | bedtools merge -i - > output/${PREFIX}.~{suffix}.mapped.merged.bed
        mv output/contig.bed output/${PREFIX}.~{suffix}.mapped.bed
        bedtools sort -i output/ref.bed | bedtools merge -i - > output/${PREFIX}.~{suffix}.ref.merged.bed
        mv output/ref.bed output/${PREFIX}.~{suffix}.ref.bed
    >>> 
    runtime {
        docker: dockerImage
        memory: memSize + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSize + " SSD"
        preemptible : preemptible
    }
    output {
        File mappedBlocksBed = glob("output/*.${suffix}.mapped.bed")[0]
        File mergedMappedBlocksBed = glob("output/*.${suffix}.mapped.merged.bed")[0]
        #File mappedLowMQBlocksBed = glob("output/*.${suffix}.mq_lq20.bed")[0]
        File unmappedBlocksBed = glob("output/*.${suffix}.unmapped.bed")[0]
        #File skippedBlocksBed = glob("output/*.${suffix}.skipped.bed")[0]
        File refBlocksBed = glob("output/*.${suffix}.ref.bed")[0]
        File mergedRefBlocksBed = glob("output/*.${suffix}.ref.merged.bed")[0]
    }
}

