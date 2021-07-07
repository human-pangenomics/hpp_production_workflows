#!/bin/bash

prefix=$1
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

## Find regions supported by each of the platforms
for platform in "hifi" "ont";do
	bedtools subtract -a $prefix.asm.bed -b $prefix.${platform}.bed | bedtools merge -d 100 -i - > $prefix.$platform.low_high.bed
	bedtools subtract -a $prefix.$platform.low_high.bed -b $prefix.asm.ends.bed > $prefix.$platform.low_high.trim1k.bed
	bedtools subtract -a $prefix.asm.bed -b $prefix.$platform.low_high.bed > $prefix.$platform.support.bed
done

## Find the blocks that are supported by at least one long read platform (ONT or HiFi)
acc $prefix.gaps.bed $prefix.{hifi,ont}.support.bed > $prefix.acc.bed 2> $prefix.acc.log
cat $prefix.acc.bed | awk '$4>0' | bedtools sort -i - | bedtools merge -i - > $prefix.acc.mrg.bed

## Find unreliable blocks
bedtools subtract -a $prefix.asm.bed -b $prefix.acc.mrg.bed | bedtools merge -d 100 -i - > $prefix.low_support.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.low_support.bed > $prefix.reliable.bed
bedtools subtract -a $prefix.low_support.bed -b $prefix.asm.ends.bed > $prefix.low_support.trim1k.bed

