#!/bin/bash

bed=$1
prefix=$2

fname=$(basename ${bed%.bed})

samtools faidx $prefix.fa
cat $prefix.fa.fai | awk '{print $1"\t"$2}' > $prefix.len.txt

printf "\n## Initial BED:\n"
cat $bed |
	awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# extend by 5% and merge overlaping blocks
bedtools slop -pct -b 0.05 -i $bed -g $prefix.len.txt | bedtools sort -i - | bedtools merge -i - > ${fname}.5perc_merged.bed

printf "\n## Extend by 5%% and merge overlaping blocks:\n"
cat ${fname}.5perc_merged.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# remove blocks shorter than 100
cat ${fname}.5perc_merged.bed | awk '{if ($3-$2 > 100) print $0}' > ${fname}.5perc_merged.gt100.bed

printf "\n## Remove blocks shorter than 100:\n"
cat ${fname}.5perc_merged.gt100.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# include whole contig if there is a block spanning over more than 50% of the contig
awk -v prefix="$prefix" -v fname="$fname" 'FILENAME == prefix".fa.fai" {len[$1]=$2} FILENAME == fname".5perc_merged.gt100.bed" {if (($3 - $2) > 0.5 * len[$1]) print $1"\t"1"\t"len[$1]; else print $0}' $prefix.fa.fai ${fname}.5perc_merged.gt100.bed | bedtools sort -i - | bedtools merge -i - > ${fname}.5perc_merged.gt100.50perc_contig.bed

printf "\n## Include the whole contig if there is a block spanning over more than 50%% of the contig:\n"
cat ${fname}.5perc_merged.gt100.50perc_contig.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# extend edge blocks (closer than 10% to one of the contig ends) by 100%
awk -v prefix="$prefix" -v fname="$fname" 'FILENAME == prefix".fa.fai" {len[$1]=$2} FILENAME == fname".5perc_merged.gt100.50perc_contig.bed" {if (($3 + ($3 - $2) * 0.1) > len[$1] || $3 < ($3 - $2) * 0.1) print $0}' $prefix.fa.fai ${fname}.5perc_merged.gt100.50perc_contig.bed | bedtools slop -pct -b 1 -i - -g $prefix.len.txt | bedtools sort -i - | bedtools merge -i - >  ${fname}.5perc_merged.gt100.50perc_contig.just_edge_blocks.bed

printf "\n## Extend edge blocks (closer than 10%% to one of the contig ends) by 100%% (Just edge blocks!)\n"
cat ${fname}.5perc_merged.gt100.50perc_contig.just_edge_blocks.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# extend middle blocks by 50% on each end
awk -v prefix="$prefix" -v fname="$fname" 'FILENAME == prefix".fa.fai" {len[$1]=$2} FILENAME == fname".5perc_merged.gt100.50perc_contig.bed" {if (($3 + ($3 - $2) * 0.1) <= len[$1] || $3 >= ($3 - $2) * 0.1) print $0}' $prefix.fa.fai ${fname}.5perc_merged.gt100.50perc_contig.bed | bedtools slop -pct -b 0.5 -i - -g $prefix.len.txt | bedtools sort -i - | bedtools merge -i - >  ${fname}.5perc_merged.gt100.50perc_contig.just_middle_blocks.bed

printf "\n## Extend middle blocks by 50%% on each end: (Just middle blocks!)\n"
cat ${fname}.5perc_merged.gt100.50perc_contig.just_middle_blocks.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

# merge all
cat ${fname}.5perc_merged.gt100.50perc_contig.just_edge_blocks.bed ${fname}.5perc_merged.gt100.50perc_contig.just_middle_blocks.bed | bedtools sort -i - | bedtools merge -i - > ${fname}.5perc_merged.gt100.50perc_contig.extended.bed

printf "\n## Merge all:\n"
cat ${fname}.5perc_merged.gt100.50perc_contig.extended.bed |
        awk '{sum += $3 -$2} END {printf "\tnumber/length of blocks: %d/ %.2f Mb\n", NR, sum/1e6}'

