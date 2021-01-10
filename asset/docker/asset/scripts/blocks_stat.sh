#!/bin/bash

prefix=$1
#centromere_bed=$2

printf "\t%45s%30s%30s\n" "HiC" "HiFi" "ONT"
printf "## Unreliable blocks:\n"
printf "\t%-20s" "Total Length ( Number )"
for platform in "hic" "hifi" "ont" "bionano"; do 
	cat $prefix.$platform.low_high.trim1k.bed | \
		awk -F '\t' '{len_blocks += $3 - $2} END {printf "%20.2f ( %d )", len_blocks/1e6, NR}'
done

printf "\n"

# intersect with centromeric regions and print the stats
printf "\n\t%45s%30s%30s\n" "HiC" "HiFi" "ONT"
printf "## Contribution of each platform\n"
printf "\t%-20s" "Total Length ( Number )"
for platform in "hic" "hifi" "ont" "bionano"; do
        bedtools intersect -a $prefix.$platform.low_high.trim1k.bed -b $prefix.low_support.trim1k.bed | \
                awk -F '\t' '{len_blocks += $3 - $2} END {printf "%20.2f ( %d )", len_blocks/1e6, NR}'
done
printf "\n"

HIC=`cat $prefix.hic.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

HIFI=`cat $prefix.hifi.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

ONT=`cat $prefix.ont.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`


mkdir -p intersect

bedtools intersect -a $prefix.hic.low_high.trim1k.bed -b $prefix.hifi.low_high.trim1k.bed > intersect/$prefix.hic_and_hifi.low_high.trim1k.bed
bedtools intersect -a $prefix.hifi.low_high.trim1k.bed -b $prefix.ont.low_high.trim1k.bed > intersect/$prefix.hifi_and_ont.low_high.trim1k.bed
bedtools intersect -a $prefix.ont.low_high.trim1k.bed -b $prefix.hic.low_high.trim1k.bed  > intersect/$prefix.ont_and_hic.low_high.trim1k.bed
bedtools intersect -a $prefix.hic.low_high.trim1k.bed \
	-b intersect/$prefix.hifi_and_ont.low_high.trim1k.bed  > intersect/$prefix.hic_and_hifi_and_ont.low_high.trim1k.bed


HIC_HIFI=`cat intersect/$prefix.hic_and_hifi.low_high.trim1k.bed | 
	awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

HIFI_ONT=`cat intersect/$prefix.hifi_and_ont.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

ONT_HIC=`cat intersect/$prefix.ont_and_hic.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

HIC_HIFI_ONT=`cat intersect/$prefix.hic_and_hifi_and_ont.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

awk -v HIC="$HIC" -v ONT="$ONT" -v HIFI="$HIFI" \
	-v HIC_HIFI="$HIC_HIFI" -v HIFI_ONT="$HIFI_ONT" -v ONT_HIC="$ONT_HIC" \
	-v HIC_HIFI_ONT="$HIC_HIFI_ONT" \
	'BEGIN{printf "\n## For Venn diagram \n\t\t(HiC,\t HiFi,\t HiC_HiFi,\t ONT,\t ONT_HiC,\t ONT_HiFi,\t HiC_HiFi_ONT): \n\t\t(%.2f,\t %.2f,\t %.2f,\t\t %.2f,\t %.2f,\t\t %.2f,\t\t %.2f)\n", \
	HIC - HIC_HIFI - HIC_ONT + HIC_HIFI_ONT, \
	HIFI - HIC_HIFI - HIFI_ONT + HIC_HIFI_ONT, \
	HIC_HIFI - HIC_HIFI_ONT, \
	ONT - ONT_HIC - HIFI_ONT + HIC_HIFI_ONT, \
	ONT_HIC - HIC_HIFI_ONT, \
	HIFI_ONT - HIC_HIFI_ONT, \
	HIC_HIFI_ONT}'

#BioNano Intersections


BN=`cat $prefix.bionano.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

mkdir -p intersect/bionano

bedtools intersect -a $prefix.bionano.low_high.trim1k.bed -b $prefix.hic.low_high.trim1k.bed  > intersect/bionano/$prefix.bionano_and_hic.low_high.trim1k.bed
bedtools intersect -a $prefix.bionano.low_high.trim1k.bed -b $prefix.hifi.low_high.trim1k.bed > intersect/bionano/$prefix.bionano_and_hifi.low_high.trim1k.bed
bedtools intersect -a $prefix.bionano.low_high.trim1k.bed -b $prefix.ont.low_high.trim1k.bed  > intersect/bionano/$prefix.bionano_and_ont.low_high.trim1k.bed

bedtools intersect -a $prefix.bionano.low_high.trim1k.bed \
        -b intersect/$prefix.hic_and_hifi.low_high.trim1k.bed  > intersect/bionano/$prefix.bionano_and_hic_and_hifi.low_high.trim1k.bed
bedtools intersect -a $prefix.bionano.low_high.trim1k.bed \
        -b intersect/$prefix.hifi_and_ont.low_high.trim1k.bed  > intersect/bionano/$prefix.bionano_and_hifi_and_ont.low_high.trim1k.bed
bedtools intersect -a $prefix.bionano.low_high.trim1k.bed \
        -b intersect/$prefix.ont_and_hic.low_high.trim1k.bed   > intersect/bionano/$prefix.bionano_and_ont_and_hic.low_high.trim1k.bed

bedtools intersect -a $prefix.bionano.low_high.trim1k.bed \
        -b intersect/$prefix.hic_and_hifi_and_ont.low_high.trim1k.bed  > intersect/bionano/$prefix.all.low_high.trim1k.bed

BN_HIC=`cat intersect/bionano/$prefix.bionano_and_hic.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

BN_HIFI=`cat intersect/bionano/$prefix.bionano_and_hifi.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

BN_ONT=`cat intersect/bionano/$prefix.bionano_and_ont.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`


BN_HIC_HIFI=`cat intersect/bionano/$prefix.bionano_and_hic_and_hifi.low_high.trim1k.bed | \
	awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`
BN_HIFI_ONT=`cat intersect/bionano/$prefix.bionano_and_hifi_and_ont.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`
BN_ONT_HIC=`cat intersect/bionano/$prefix.bionano_and_ont_and_hic.low_high.trim1k.bed | \
        awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`

ALL=`cat intersect/bionano/$prefix.all.low_high.trim1k.bed | \
	awk -F '\t' '{len_blocks += $3 - $2} END {printf "%.4f", len_blocks/1e6}'`


awk -v BN_HIC="$BN_HIC" -v BN_HIFI="$BN_HIFI" -v BN_ONT="$BN_ONT" \
	-v BN_HIC_HIFI="$BN_HIC_HIFI" -v BN_HIFI_ONT="$BN_HIFI_ONT" -v BN_ONT_HIC="$BN_ONT_HIC" \
	-v ALL="$ALL" -v BN="$BN"\
        'BEGIN{printf "\n## BIONANO\n\t\t(HiC,\t HiFi,\t HiC_HiFi,\t ONT,\t ONT_HiC,\t ONT_HiFi,\t ALL,\t None): \n\t\t(%.2f,\t %.2f,\t %.2f,\t\t %.2f,\t %.2f,\t\t %.2f,\t\t %.2f,\t %.2f)\n", \
        BN_HIC - BN_HIC_HIFI - BN_HIC_ONT + ALL, \
        BN_HIFI - BN_HIC_HIFI - BN_HIFI_ONT + ALL, \
        BN_HIC_HIFI - ALL, \
        BN_ONT - BN_ONT_HIC - BN_HIFI_ONT + ALL, \
        BN_ONT_HIC - ALL, \
        BN_HIFI_ONT - ALL, \
        ALL, \
	BN - (BN_HIC + BN_HIFI + BN_ONT - BN_HIC_HIFI - BN_HIFI_ONT - BN_ONT_HIC + 2 * ALL)}'


