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

for platform in "hic" "hifi" "ont" "bionano";do
	bedtools subtract -a $prefix.asm.bed -b $prefix.${platform}.bed | bedtools merge -d 100 -i - > $prefix.$platform.low_high.bed
	bedtools subtract -a $prefix.$platform.low_high.bed -b $prefix.asm.ends.bed > $prefix.$platform.low_high.trim1k.bed
	bedtools subtract -a $prefix.asm.bed -b $prefix.$platform.low_high.bed > $prefix.$platform.support.bed
done
acc $prefix.gaps.bed $prefix.{hic,hifi,ont}.support.bed > $prefix.acc.bed 2> $prefix.acc.log
awk '$4>1' $prefix.acc.bed | bedtools merge -i - > $prefix.acc.gt2.mrg.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.acc.gt2.mrg.bed | bedtools merge -d 100 -i - > $prefix.low_support.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.low_support.bed > $prefix.reliable.bed
bedtools subtract -a $prefix.low_support.bed -b $prefix.asm.ends.bed > $prefix.low_support.trim1k.bed

awk '$4>2' $prefix.acc.bed | bedtools merge -i - > $prefix.acc.gt3.mrg.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.acc.gt3.mrg.bed | bedtools merge -d 100 -i - > $prefix.gt3.low_support.bed
bedtools subtract -a $prefix.asm.bed -b $prefix.gt3.low_support.bed > $prefix.gt3.reliable.bed
bedtools subtract -a $prefix.gt3.low_support.bed -b $prefix.asm.ends.bed > $prefix.gt3.low_support.trim1k.bed
