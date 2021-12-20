while getopts a:b:o:t:m:h flag
do
    case "${flag}" in
	    h) printf "\nThis command is for producing the final bed output of the CopyVal pipeline\nFor more information check https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage\n\tRun command:\n\tbash combine_alt_filtered_beds.sh \ \n\t -a <BEDS_TAR_GZ> \ \n\t -b <ALT_FILTERED_BEDS_TAR_GZ> \ \n\t -m <COLORS_TEXT> \ \n\t -t <TRACK_NAME> \ \n\t -o <OUTPUT_BED>\n"; exit 1 ;;
        a) correctedBedsTarGz=${OPTARG};;
        b) altRemovedBedsTarGz=${OPTARG};;
        o) outputBed=${OPTARG};;
	t) trackName=${OPTARG};;
	m) mapColorsText=${OPTARG};;
    esac
done

mkdir corrected altRemoved
tar --strip-components 1 -xvzf ${correctedBedsTarGz} --directory corrected
tar --strip-components 1 -xvzf ${altRemovedBedsTarGz} --directory altRemoved
         
ERROR_BED_1=corrected/$(ls corrected | grep "error")
DUPLICATED_BED_1=corrected/$(ls corrected | grep "duplicated")
HAPLOID_BED_1=corrected/$(ls corrected | grep "haploid")
COLLAPSED_BED_1=corrected/$(ls corrected | grep "collapsed")
 
ERROR_BED_2=altRemoved/$(ls altRemoved | grep "error")
DUPLICATED_BED_2=altRemoved/$(ls altRemoved | grep "duplicated")
HAPLOID_BED_2=altRemoved/$(ls altRemoved | grep "haploid")
COLLAPSED_BED_2=altRemoved/$(ls altRemoved | grep "collapsed")

mkdir beds
## Haploid to Error
bedtools intersect \
    -a ${HAPLOID_BED_1} \
    -b ${ERROR_BED_2}  \
    > beds/Eh.bed
## Haploid to Duplicated
bedtools intersect \
    -a ${HAPLOID_BED_1} \
    -b ${DUPLICATED_BED_2} \
    > beds/Dh.bed
## Collapsed to Haploid
bedtools intersect \
    -a ${COLLAPSED_BED_1} \
    -b ${HAPLOID_BED_2} \
    > beds/Hc.bed
## Collapsed to Error
bedtools intersect \
    -a ${COLLAPSED_BED_1} \
    -b ${ERROR_BED_2} \
    > beds/Ec.bed
## All regions with their components changed
cat beds/Eh.bed beds/Dh.bed beds/Hc.bed beds/Ec.bed | \
    bedtools sort -i - | \
    bedtools merge -i - \
    > beds/changed.bed
## Haploid
bedtools subtract \
    -a ${HAPLOID_BED_1}  \
    -b beds/changed.bed \
    > beds/Hh.bed
## Collapsed
bedtools subtract \
    -a ${COLLAPSED_BED_1} \
    -b beds/changed.bed \
    > beds/Cc.bed

## Error blocks; Ee, Eh and Ec
cat ${ERROR_BED_1} | awk '{print $0"\tEe"}' >> beds/all.bed
cat beds/Eh.bed | awk '{print $0"\tEh"}' >> beds/all.bed
cat beds/Ec.bed | awk '{print $0"\tEc"}' >> beds/all.bed

## Duplicated blocks; Dd and Dh
cat ${DUPLICATED_BED_1} | awk '{print $0"\tDd"}' >> beds/all.bed
cat beds/Dh.bed | awk '{print $0"\tDh"}' >> beds/all.bed

## Haploid blocks; Hc and Hh
cat beds/Hc.bed | awk '{print $0"\tHc"}' >> beds/all.bed
cat beds/Hh.bed | awk '{print $0"\tHh"}' >> beds/all.bed

## Collapsed blocks; Cc only
cat beds/Cc.bed | awk '{print $0"\tCc"}' >> beds/all.bed


echo "track name=\"${trackName}\" visibility=1 itemRgb="On"" > ${outputBed}
bedtools sort -i beds/all.bed > beds/all.sorted.bed
awk 'FNR==NR{c[$1]=$2;next}{print $0"\t0\t+\t"$2"\t"$3"\t"c[$4]}' \
    ${mapColorsText} \
    beds/all.sorted.bed \
    >> ${outputBed}

rm -rf altRemoved corrected beds

