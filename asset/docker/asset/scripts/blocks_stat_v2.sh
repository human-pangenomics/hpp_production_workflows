#!/bin/bash

prefix=$1
cenSatRegions=$2
printf "Unreliable blocks:\n"

printf "\n###########\nLOW SUPPORT\n###########\n"
cat $prefix.low_support.trim1k.bed | python3 ${BLOCKS_STAT_PY}

printf "\n#####################################\nLOW SUPPORT AND CENTROMERIC SATELLITE\n#####################################\n"
bedtools intersect -a $prefix.low_support.trim1k.bed -b $cenSatRegions | bedtools sort -i - | bedtools merge -i - | python3 ${BLOCKS_STAT_PY}

printf "\n#####################################\nCENTROMERIC SATELLITE\n#####################################\n"
cat $cenSatRegions | python3 ${BLOCKS_STAT_PY}

printf "\n####\nHIFI\n####\n"
cat $prefix.hifi.low_high.trim1k.bed | python3 ${BLOCKS_STAT_PY}

printf "\n###\nONT\n###\n"
cat $prefix.ont.low_high.trim1k.bed | python3 ${BLOCKS_STAT_PY}
