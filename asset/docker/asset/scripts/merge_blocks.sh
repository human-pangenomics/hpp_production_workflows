#!/bin/bash

prefix=$1
cenSatRegions=$2
# Get gaps
detgaps $prefix.fa > $prefix.gaps.bed
 
# Get scaffold length
samtools faidx $prefix.fa
 
# Get scaffold region
awk '{print $1"\t0\t"$2}' $prefix.fa.fai > $prefix.asm.bed
 
# Get scaffold ends while ignoring scaffolds <2kb
cat $prefix.asm.bed | \
 awk '($3-$2) > 2000 {print $1"\t0\t1000\n"$1"\t"($3-1000)"\t"$3}' \
 > $prefix.asm.ends.bed

## Find regions supported by each of 4 platforms
for platform in "hifi" "ont" "bionano";do
	bedtools subtract -a $prefix.asm.bed -b $prefix.${platform}.bed | bedtools merge -d 100 -i - > $prefix.$platform.low_high.bed 
	cat $prefix.$platform.low_high.bed | awk '{if (($3 -$2) > 10) { print $0 }}' > $prefix.$platform.low_high.gt10.bed
	bedtools subtract -a $prefix.$platform.low_high.bed -b $prefix.asm.ends.bed > $prefix.$platform.low_high.trim1k.bed
	cat $prefix.$platform.low_high.trim1k.bed | awk '{if (($3 -$2) > 10) { print $0 }}' > $prefix.$platform.low_high.trim1k.gt10.bed
	bedtools subtract -a $prefix.asm.bed -b $prefix.$platform.low_high.gt10.bed > $prefix.$platform.support.bed
done

## Since HiC reads and Bionano molecules are hardly aligned to centromeric regions, these regions are added to the regions supported by HiC and Bionano
for platform in "bionano";do
	bedtools subtract -a $prefix.$platform.low_high.bed -b $cenSatRegions > $prefix.$platform.low_high.no_censat.bed
	cat $prefix.$platform.low_high.no_censat.bed | awk '{if (($3 -$2) > 10) { print $0 }}' > $prefix.$platform.low_high.no_censat.gt10.bed
	bedtools subtract -a $prefix.$platform.low_high.trim1k.bed -b $cenSatRegions > $prefix.$platform.low_high.trim1k.no_censat.bed
	cat $prefix.$platform.low_high.trim1k.no_censat.bed | awk '{if (($3 -$2) > 10) { print $0 }}' > $prefix.$platform.low_high.trim1k.no_censat.gt10.bed
        bedtools subtract -a $prefix.asm.bed -b $prefix.$platform.low_high.no_censat.gt10.bed > $prefix.$platform.support.with_censat.bed
done

## Find the blocks that are supported by at least one long read platform (ONT or HiFi)
#acc $prefix.gaps.bed $prefix.{hifi,ont}.support.bed $prefix.hic.support.with_censat.bed > $prefix.acc_1.bed 2> $prefix.acc_1.log
acc $prefix.gaps.bed $prefix.{hifi,ont}.support.bed $prefix.bionano.support.with_censat.bed > $prefix.acc_2.bed 2> $prefix.acc_2.log
cat $prefix.acc_2.bed | awk '$4>1' | bedtools sort -i - | bedtools merge -i - > $prefix.acc.gt2.mrg.bed

## Find unreliable blocks
bedtools subtract -a $prefix.asm.bed -b $prefix.acc.gt2.mrg.bed | bedtools merge -d 100 -i - > $prefix.low_support.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.low_support.bed > $prefix.reliable.bed
bedtools subtract -a $prefix.low_support.bed -b $prefix.asm.ends.bed > $prefix.low_support.trim1k.bed
cat $prefix.low_support.trim1k.bed | awk '{if (($3 -$2) > 10) { print $0 }}' > $prefix.low_support.trim1k.gt10.bed

