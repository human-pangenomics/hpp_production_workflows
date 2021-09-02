## Heterozygosity of HPRC samples

This page describes the pipeline that produces the neccessary statistics for estimating the heterozygosities of the HPRC-Year1 samples.

### Dipcall outputs are the main inputs:

As a part of the assembly QC pipeline, Dipcall is used for calling the variants in each HPRC sample. Dipcall is a reference-based variant caller 
and GRCh38 has been used as the reference for the current results. 
The heterozygosity pipeline needs the Dipcall outputs and these are four the main files neccessary for running the current pipeline:

1. `${SAMPLE}.f1_assembly_v2_genbank.dip.vcf.gz`: The VCF file that contains the called variants.
2. `${SAMPLE}.f1_assembly_v2_genbank.hap1.bam`: The BAM file that contains the alignment of the paternal assembly to the reference.
3. `${SAMPLE}.f1_assembly_v2_genbank.hap2.bam`: The BAM file that contains the alignment of the maternal assembly to the reference.
4. `${SAMPLE}.f1_assembly_v2_genbank.dip.bed`: The BED file that points to the confident regions in the reference for variant calling. More information about the selection of confident regions is available [here](https://github.com/lh3/dipcall).

### Parsing the variants and calculating the main statistics
The main statistics are calculated by the passing the VCF file to the script `het_estimator.py`. This script is available in `docker/scripts`.
````
# Filter the variants not in the confident regions
bedtools intersect -header -a <(zcat ${vcf_gz}) -b ${confident_bed} > confident.vcf

# Produce the statistics
python3 het_estimator.py --vcf confident.vcf > het_stats.txt
````

`het_stats.txt` is a tab-delimited file like below:
```
#region	het_snp	hom_non_ref_snp	het_snp_within_insertion	insertion_hap1	insertion_hap2
autosome	2289811	1597847	7155	6830105	6448089
par	3848	2582	98	60063	34494
non_par_X	0	0	0	0	0

Het Ratio (het_snp / hom_non_ref_snp) in autosome and par: 1.43
```

It contains 7 columns:

1. `region`: It can be one of the three regions; 
    - autosomes, 
    - PAR (both chrY and chrX) and
    - non-PAR (only chrX)


     The statistics are stratified by these regions since male and female samples may have different heterozygosity patterns in the last two regions; PAR and non-PAR-chrX. 
The non-PAR of chrY is not included since we cannot find any heterozygous variants in this region for either male or female samples.

3. `het_snp`: The number of heterozygous SNPs in the confident regions.
4. `hom_non_ref_snp`: The number of non-reference homozygous SNPs in the confident regions
5. `het_snp_within_insertion`: In some cases there are insertions at the same reference location in both haplotypes. 
These insertions are with respect to the reference and the VCF file does not provide the heterozygous SNPs that may exist within these insetions. 
So to find those SNPs a global alignment is done between the inserted sequences. `edlib` python library does the global alignment and finds the path with 
the least edit distance (the scores would be; match score = 0, mismatch score = -1, gap score = -1). 

    Here is an example:

    <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/snp_within_insertion.png" width="1000" height="200">

    This is the corresponding vcf record:
    ```
    ['chr2', '120674858', '.', 'G', 'GGCGGCACGAGGGGTGCACTTTGAGCCCA,GGCAGCACGAGGGGTGCACTTTGAGCCCA', '30', '.', '.', 'GT:AD', '2|1:0,1,1']
    ```
    By doing the global alignment between the inserted sequences we will have:

    ```
    GGCAGCACGAGGGGTGCACTTTGAGCCCA
    |||.|||||||||||||||||||||||||
    GGCGGCACGAGGGGTGCACTTTGAGCCCA
    ```

    So there is one heterozygous SNP within insertions.

5. `insertion_hap1`: The number of inserted bases in the haplotype 1 (paternal) w.r.t the other haplotype. There are three cases:

    - Both haplotypes have insertions happening at the same reference location. In this case we iterate over the global alignment explained in the previous item. 
The number of inserted bases in haplotype 1 is counted and added to the column `insertion_hap1`.
        
        Here is an example:

        <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/insertion_insertion.png" width="1000" height="200">

        This is the corresponding vcf record.
        ```
        ['chr1', '144604740', '.', 'T', 'TATATATATATATATATATATATACACACACACAC,TATATATATATATATATATATATATATATACACACACACAC', '30', '.', '.', 'GT:AD', '2|1:0,1,1']
        ```

        The global alignment

        ```
        TATATATATATATATATATATATATATATACACACACACAC
        ||||||||||||||||||||||||------|||||||||||
        TATATATATATATATATATATATA------CACACACACAC
        ```
        So there is an insertion of length 6 in hap1.

    - There is a deletion in haplotype 2 and a shorter deletion (or no deletion) in haplotype 1 starting at the same reference location. In this case we simply subtract the length of deletions and add the number to the column `insertion_hap1`. 
    
        Here is an example:
        <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/deletion_deletion.png" width="1000" height="200">

    - There is an insertion in haplotype 1 and a deletion (or no deletion) in haplotype 2 starting at the same reference location. We add both the insertion and deletion length to the column `insertion_hap1`. 
        
        Here is an example:
        <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/insertion_deletion.png" width="1000" height="200">

    
6. `insertion_hap2`: This is the same as `insertion_hap1` but for haplotype 2 (maternal).

### Calculating the total length of the confident blocks in the assembly

 In order to find the total length of the assembly blocks aligned to the confident regions we can use the 
 (hap1 and hap2) BAM files and project the confident reference regions onto each of the assemblies. The python script `project_blocks.py` can do this projection:
 
 ```
 # Convert bam 2 paf by paftools.js
 k8 paftools.js sam2paf <(samtools view -h -q 5 {hap1_bam}) > hap1.paf
 
 # Project confident regions
 # hap1_confident_projection.bed contains the projection of the blocks in {confident_bed} onto the assembly (That's what we are interested in)
 # hap1_confident_projectable.bed contains the confident blocks that could be projected (reference space)
python3 project_blocks.py --mode 'ref2asm' --paf hap1.paf --blocks {confident_bed} --outputProjection hap1_confident_projection.bed --outputProjectable hap1_confident_projectable.bed 

# Sort and merge
bedtools sort -i  hap1_confident_projection.bed | bedtools merge -i - > hap1_confident_projection.merged.bed

# Get the total length
cat hap1_confident_projection.merged.bed | awk '{sum += $3 -$2} END {print sum}'
 ```

### Using the statistics to estimate the heterozygosity:

#### First approach (snp-based heterozygosity ratio)
One metric that could be found in the literature is the heterozygosity ratio [like in this paper](https://academic.oup.com/genetics/article/204/3/893/6066348), which is the ratio of the number of heterozygous snps divided by the number of non-reference homozygous snps. 
This number is being reported in the last line of `het_estimator.py` 's output. If we calculate this number for
all HPRC samples and plot them against the super population we can see the dependency.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/het_ratio_pop.png" width="600" height="400">

For this metric only the autosomal regions are used.

#### Second approach (Indel-inclusive heterozygosity percent)
The other metric worth considering is the percentage of the diploid assembly bases that are either snps or insertions. We can summarize it in the formula below:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/het_percent_formula.png" width="800" height="80">

We can again notice the dependency between the super population and this metric:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/het_estimation/het_estimation/images/het_percent_pop.png" width="600" height="400">

For this metric only the autosomal regions are used.

### Results Availability
A table is available [here](https://docs.google.com/spreadsheets/d/1SyyaY3vUzv89rCKqto9j0y5t6wE_Q8I5iMLlAJxTAbs/edit?usp=sharing) that contains the numbers and estimated values for all HPRC-Year1 samples.
