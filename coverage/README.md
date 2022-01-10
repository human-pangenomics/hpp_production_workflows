## Flagger 

### Overview
**Flagger** is a read-based pipeline that can detect different types of mis-assemblies in a draft diploid assembly. One core component of this pipeline is another pipeline named [**Partitioner**](https://github.com/human-pangenomics/hpp_production_workflows/edit/asset/coverage/docs/coverage/README.md). Partitioner recieves the read alignments to the draft assembly and partition the assembly into 4 main components; erroneous, (falsely) duplicated, haploid and collapsed.

Flagger has 7 steps:
- Align reads to the diploid assembly
- Phase the ambiguous alignments using [the phasing pipeline](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/README.md)
- Run Partitioner on the alignments
- Call variants 
- Remove the alignments with alternative alleles
- Run Paritioner using the alignments with no alternative allele
- Combine the Partitioner outputs
### 1. Align Reads
The ONT and HiFi reads can be aligned to a diploid assembly (~ 6Gbases long) with winnowmap. Since the assembly is diploid the expected base-level coverage should be half of the sequencing coverage.
Here are the main commands for producing the alignments (taken from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] -Y -L --eqx --cs -I8g <(cat pat_asm.fa mat_asm.fa) reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````
### 2. Phase reads
In this step the reads with multiple alignments are phased using single-base markers. In other words all the secondary and primary
alignments of the same read are scored based on marker consistency and 
the alignment with the highest score is selected as the primary alignment. The output of this section is 
a corrected version of the input bam file, in which the primary and secondary alignments are swapped 
whenever neccessary.

More information about the phasing pipeline is available [here](https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/docs/phasing/README.md)

### 3. Run Partitioner on the alignments
The corrected bam file is then given as input to Partitioner. Partitioner outputs a bed file for each of the 4 components; 
erroneous, duplicated, haploid and collapsed. Any component other than the haploid one is pointing to unreliable blocks in
assembly. The 4 components are explained in detail [here](https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/docs/coverage). 

More information about Partitioner is available [here](https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/docs/coverage)

### 4. Call variants 
By calling variants it is possible to detect the regions that needs polishing or the regions that have alignments from the wrong haplotype. It is recommeneded to use [Deepvariant](https://github.com/google/deepvariant) for calling variants from HiFi alignments and [Pepper-Margin-Deepvariant](https://github.com/kishwarshafin/pepper) for ONT. 
````
## For HiFi
## Taken from deepvariant doc
BIN_VERSION="1.3.0"
docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type="PACBIO" \
  --ref="/input/${ASSEMBLY_FASTA}" \
  --reads="/input/${INPUT_BAM}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --make_examples_extra_args="keep_supplementary_alignments=true, min_mapping_quality=0" \
  --call_variants_extra_args="use_openvino=true" \
  --num_shards=$(nproc) \
  --dry_run=false 
  
## For ONT
## Taken from pepper-margin-deepvariant doc
sudo docker run \
  -v ${INPUT_DIR}:/input \
  -v ${OUTPUT_DIR}:/output \
  kishwars/pepper_deepvariant:r0.6 \
  run_pepper_margin_deepvariant call_variant \
  -b "/input/${INPUT_BAM}" \
  -f "/input/${ASSEMBLY_FASTA}" \
  -o "/output" \
  -t $(nproc) \
  --ont_r9_guppy5_sup \
  --pepper_include_supplementary \
  --dv_min_mapping_quality 0 \
  --pepper_min_mapping_quality 0 \
  
 
# --ont_r9_guppy5_sup is preset for ONT R9.4.1 Guppy 5 "Sup" basecaller
# for ONT R10.4 Q20 reads: --ont_r10_q20
````

Note that for both variant callers, the minimum mapping quality is set to 0 which is neccessary to do if the assembly under evaluation is diploid.

### 5. Remove the alignments with alternative alleles
The called variants are then filtered to include only the biallelic snps.
````
## Get the biallelic snps
bcftools view -Ov -f PASS -m2 -M2 -v snps ${OUTPUT_VCF} > ${SNPS_VCF}
````
By having the biallelic snps it is possible to find the alignments with alternative alleles, remove them from the bam file and produce a new bam file.
`filter_alt_reads` is a program that can be used for this aim.
```
## Run filter_alt_reads to get a bam file with no alternative-contained alignments
docker run \
 -v ${INPUT_DIR}:/input \
 -v ${OUTPUT_DIR}:/output \
 quay.io/masri2019/hpp_coverage:latest \
 filter_alt_reads \
 -i "/input/${INPUT_BAM}" \
 -o "/output/${ALT_FILTERED_BAM}"
 -f "/output/${ALT_BAM}"
 -v "${SNPS_VCF}"
 -t $(nproc)
```
`${ALT_FILTERED_BAM}` is the bam file with no alignments that contain alternative alleles and `${ALT_BAM}` includes the removed alignments.

### 6. Run Paritioner on the alignments with no alternative allele
`${ALT_FILTERED_BAM}` produced in the previous step will be used as a new input for Partitioner. This step is same as the 3rd step except that the input bam file here does not have the alternative-contained alignments.


### 7. Combine the Partitioner outputs in steps 3 and 6
The partitioner outputs in steps 3 and 6 are expected to have a huge overlap but they are not the same. As an example one region that was detected as collapsed in step 3 may be categorized as haploid in step 6. This component change is showing that the flagged region is assembled correctly but 
its homologous region in the other haplotype is not assembled correctly and that's why the reads from the other haplotype aligned there. 
So by combining these two sets of partitioner output it is possible to infer more information about the unreliable blocks in the assembly.

The patitioner output is a gzipped tar file that contain 4 bed files one for each component.
```
bash combine_alt_removed_beds \
-a ${BEDS_TAR_GZ_STEP3} \
-b ${BEDS_TAR_GZ_STEP6} \
-m /home/scripts/colors.txt \
-t ${SAMPLE_NAME} \
-o ${OUTPUT_BED}
```

In `${OUTPUT_BED}` one of the 8 components below is assigned to each block. 

### Components
Only `Hh` and `Hc` point to the regions with expected read support. `Hc` also shows a mis-assembly in another haplotype.
|Component|Initial|After|Color |Description|
|:--------|:------|:----|:-----|:----------|
|Cc |Collapsed |**Collapsed** |Purple| Two highly similar haplotypes are collapsed into this block |
|Hc  |Collapsed | **Haploid** |Blue|This block is assembled correctly. It also has false alignments from a not assembled haplotype |
|Dd  |Duplicated |**Duplicated** |Orange| This block is a false duplication of another block. (Mainly low-MAPQ alignments with half of the expected coverage)|
|Ee  |Erroneous |**Erroneous** |Dark Red| This block has low read coverage |
|Dh  |Haploid |**Duplicated** |Yellow| This block is a false duplication of another block like `Dd`, it also has false alignments from a not assembled haplotype |
|Eh  |Haploid| **Erroneous** |Red| This block needs polishing |
|Hh  |Haploid| **Haploid** | Green|This block is correctly assembled and has the expected read coverage |
|Ec  |Collapsed| **Errorneous** | Pink|This block needs polishing. It also has alignments from multiple not-assembled haplotypes and after removing the false alignments it does not have the expected read coverage|

`Initial` column shows the component the block has been assigned to before removing the alignments with alternative alleles and `After` shows the component after removing. Each of these components has their own color when they are shown in the IGV or Genome Browser.


### Data, Source Code and Workflows Availability

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

We have used the Genbank version of the HPRC-Y1 assemblies.

The Python scripts, C source codes and the binary files are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/docker/coverage/scripts

The wdl files that have been used for this analysis are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/wdl/tasks


### Results Availability
### ATTENTION: These Results ARE OLD! Will be updated soon.

The results are available in 

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/

The most updated results are under the `V1` directory. It includes the results of the HiFi-based coverage analysis.
The directory structure for the results of each sample is as below:

````
HG002
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.cov_dist.pdf
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.table
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.tables.tar.gz
├── combined
├── contig_based
├── filtered
├── dup_corrected
└── whole_genome_based
````

The pdf file contains the plots of the coverage distributions for the whole assembly and the contigs longer than 5Mb which has been used for the contig-based analysis. Different model components are also drawn with different colors.

The table files contain the models fit to the distributions. These models are used for finding the thresholds and categorizing the regions into one of the four components; error, duplicated, haploid and collapsed.

The four subdirectories are having the categorized blocks. They are the results of different parts of the pipeline and named based on this.
For example `combined` directory contains the combination of `contig_based` and `whole_genome_based` blocks (section 6) and `dup_corrected` contains the corrections explained in section 7.
Each directory has four bed files and they are named based on the component they are pointing to.

````
filtered/
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.collapsed.bed
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.duplicated.bed
├── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.error.bed
└── HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.filtered.haploid.bed
````

For example the bed file that ends with `.haploid.bed` is pointing to the haploid blocks.

It is recommended to use the bed files in the `filtered` subdirectory.


