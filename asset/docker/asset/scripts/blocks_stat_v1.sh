#!/bin/bash

prefix=$1
cenSatRegions=$2
printf "Unreliable blocks:\n"

printf "\n###########\nLOW SUPPORT\n###########\n"
cat $prefix.low_support.trim1k.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n#####################################\nLOW SUPPORT AND CENTROMERIC SATELLITE\n#####################################\n"
bedtools intersect -a $prefix.low_support.trim1k.gt10.bed -b $cenSatRegions | python3 ${BLOCK_STATS_PY}

printf "\n####\nHIFI\n####\n"
cat $prefix.hifi.low_high.trim1k.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n###\nONT\n###\n"
cat $prefix.ont.low_high.trim1k.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n###\nHIC\n###\n"
cat $prefix.hic.low_high.trim1k.no_censat.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n#######\nBIONANO\n#######\n"
cat $prefix.bionano.low_high.trim1k.no_censat.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n############\nHIFI AND ONT\n############\n"
bedtools intersect -a $prefix.hifi.low_high.trim1k.gt10.bed -b $prefix.ont.low_high.trim1k.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n###############\nHIC AND BIONANO\n###############\n"
bedtools intersect -a $prefix.hic.low_high.trim1k.no_censat.gt10.bed -b $prefix.bionano.low_high.trim1k.no_censat.gt10.bed > $prefix.hic_and_bionano.low_high.trim1k.no_censat.gt10.bed
cat $prefix.hic_and_bionano.low_high.trim1k.no_censat.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n########################\nHIFI AND HIC AND BIONANO\n########################\n"
bedtools intersect -a $prefix.hifi.low_high.trim1k.gt10.bed -b $prefix.hic_and_bionano.low_high.trim1k.no_censat.gt10.bed | python3 ${BLOCK_STATS_PY}

printf "\n#######################\nONT AND HIC AND BIONANO\n#######################\n"
bedtools intersect -a $prefix.ont.low_high.trim1k.gt10.bed -b $prefix.hic_and_bionano.low_high.trim1k.no_censat.gt10.bed | python3 ${BLOCK_STATS_PY}

