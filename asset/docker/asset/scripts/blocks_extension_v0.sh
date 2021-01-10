#!/bin/bash

bed=$1
prefix=$2

fname=$(basename ${bed%.bed})

samtools faidx $prefix.fa
cat $prefix.fa.fai | awk '{print $1"\t"$2}' > $prefix.len.txt

printf "\n## Initial BED:\n"
cat $bed |
	awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# extend all blocks by 50% on each end
cat $bed | bedtools slop -pct -b 0.5 -i - -g $prefix.len.txt | bedtools sort -i - | bedtools merge -i - >  ${fname}.50per.extended.bed

printf "\n## Extend 50%%:\n"
cat ${fname}.50per.extended.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'


