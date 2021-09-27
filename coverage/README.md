## A Mixture Model-Based Coverage Analysis For HPRC Y1 Assemblies

### Overview
On this page we are explaining the steps toward the coverage analysis of HPRC-Y1 assemblies. This analysis is still under development and this page will be updated once there is a newer version of the results. Here is a summary of this analysis 
- Align the long reads to each assembly
- Calculate the read coverage of each base of the assembly 
- Fit a mixture model to the coverage distribution
- Extract the blocks assigned to the model's 4 main components: erroreous, duplicated, haploid, and collapsed.

### 1. Read Alignment
The ONT and HiFi reads are aligned to each diploid assembly (~ 6Gbases long) with winnowmap v2.03. Since we are aligning the reads to the diploid assembly the expected base-level coverage should be half of the sequencing coverage.
Here are the main commands used for producing the alignments (taken from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] -Y -L --eqx --cs -I8g <(cat pat_asm.fa mat_asm.fa) reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````


### 2. Calculating Depth of Coverage
With `samtools depth` we could calculate the the depth of coverage for each assembly base. The output of `samtools depth -aa` is like below. (`-aa` 
option allows reporting the bases with no coverage)
````
contig_1  1 0
contig_1  2 1
contig_1  3 1
contig_1  4 1
contig_1  5 2
contig_1  6 2
````
In the order in which they appear, the columns are showing the contig name, the base coordinate and the coverage. 
Each base has a separate line even if consecutive bases are having the same coverage. 
In order to make this output more efficient we coverted it to the format below
````
>contig_1 6
1 1 0
2 4 1
5 6 2
````
Where each contig's name appears only once before the first block of that contig and the consecutive bases with the same coverage take only one line.
The number that comes after the name of the contig is the contig size. We added the suffix `.cov` to the files with this format.
Here are the main commands for producing the coverage files:
````
samtools depth -aa -Q 0 read_alignment.bam > read_alignment.depth
./depth2cov -d read_alignment.depth -f asm.fa.fai -o read_alignment.cov
````
`depth2cov` is a binary executable that converts the output of `samtools depth` to `.cov` format. Its source code is written in C and is available in `docker/coverage/scripts`

Since the reads aligned to the homozygous regions are expected to have low mapping qualities we don't filter reads based on their mapping qualities.
In the figure below you can see the histograms of mapping qualities and the distributions of alignment indentities for HG00438 as an example. Three sets of alignments are shown here; the alignments to the diploid assembly and to each haploid assembly (maternal and paternal) separately. The alignments to the haploid assemblies are shown here just for comparison and as it was mentioned above they are not used for the current analysis.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_mapq_hist.png" width="700" height="275">

About 20% of the diploid alignments are having MAPQs lower than 20. One of the future developments will be on accurately phasing the reads with low MAPQ alignments. 

### 3. Coverage Distribution and Fitting The Model


The frequencies of coverages can be calculated w/ `cov2counts` which is a binary executable. The source code is written in C and is available in `docker/coverage/scripts`
````
./cov2counts -i read_alignment.cov -o read_alignment.counts
````
The output file with `.counts` suffix is a 2-column tab-delimited file; the first column shows coverages and the second column shows the frequencies of
those coverages. Using `.counts` files we can produce distribution plots easily. For example below we are showing the coverage distribution for HG00438 diploid assembly.
<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_dip_hifi_cov_dist.png" width="700" height="275">

Then we pass each `.counts` file to `fit_model_extra.py` whose job is to fit the mixture model and find the best parameters through 
[Expectation Maximization (EM)](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm). This mixture model
consists of 4 main components and each component represents a specific type of regions.

The 4 components:

1. **Error component**, which is modeled by a poisson distribution. To avoid overfitting, this mode only uses the coverages below 10 so its mean is limited to be between 0 and 10. It represents the regions with very low read support.
2. **Duplicated component**, which is modeled by a gaussian distribution whose mean is constrained to be half of the haploid component's mean. It should mainly represents the falsely duplicated regions and the weight of this component is usually near to zero. It is worth noting that according to the recent 
[T2T paper, The complete sequence of a human genome,](https://www.biorxiv.org/content/10.1101/2021.05.26.445798v1.abstract) there exist
 some satellite arrays (especially HSAT1) where the ONT and HiFi coverage drops systematically due to bias in sample preparation and sequencing.
 As a result this mode should contain a mix of duplicated and coverage-biased blocks, which are not easy to be distiguished through the current analysis.
3. **Haploid component**, which is modeled by a gaussian distribution. It represents blocks with the coverages that we expect for the blocks of an error-free assembly.
4. **Collpased component**, which is actually a set of components each of which follows a gaussian distribution and their means are constrained to be
multiples of the haploid component's mean.

Here is the command that fits the model:
````
python3 fit_model_extra.py --counts read_alignment.counts --output read_alignment.table
````

Its output is a file with `.table` suffix. It contains a TAB-delimited table with the 7 fields described below:
|Col|Type  |Description                               |
|:--|:----:|:-----------------------------------------|
|coverage  |int|A depth of coverage                       |
|freq  |float   |The frequency of the coverage value in the first column                     |
|fit  |float   |The frequency value fit to the model    |
|error  |float   |The weight of the error component       |
|duplicated  |float  |The weight of the duplicated component               |
|haploid  |float|The weight of the haploid component                      |
|collapsed  |float   |The weight of the collapsed component                    |

Here is an example of such a table:
````
#coverage	freq	fit	error	duplicated	haploid	collapsed
0	1633030	2479818.7542	0.9973	0.0009	0.0019	0.0000
1	552306	1014768.3370	0.9890	0.0044	0.0067	0.0000
2	311274	425804.6727	0.9564	0.0207	0.0229	0.0000
3	241294	195958.8474	0.8434	0.0856	0.0710	0.0000
4	201587	117491.9394	0.5708	0.2616	0.1676	0.0000
5	207799	108998.0736	0.2497	0.4973	0.2530	0.0000
6	214003	141291.2687	0.0782	0.6512	0.2706	0.0000
7	280495	207289.2292	0.0216	0.7251	0.2533	0.0000
8	317883	309510.6691	0.0059	0.7635	0.2306	0.0000
...
...
...
60	1547886	1045679.3816	0.0000	0.0000	0.7351	0.2649
61	1418178	920434.1445	0.0000	0.0000	0.6731	0.3269
62	1326917	819701.9290	0.0000	0.0000	0.6031	0.3969
63	1285300	740624.2274	0.0000	0.0000	0.5275	0.4725
64	1211029	680304.0405	0.0000	0.0000	0.4493	0.5507
65	1111166	635906.0944	0.0000	0.0000	0.3724	0.6276
66	1034004	604733.6596	0.0000	0.0000	0.3004	0.6996
67	948793	584283.3204	0.0000	0.0000	0.2362	0.7638
68	849748	572280.0759	0.0000	0.0000	0.1814	0.8186
69	800708	566695.7729	0.0000	0.0000	0.1364	0.8636
70	766435	565754.1385	0.0000	0.0000	0.1008	0.8992
````

In the figure below, we are showing the model components. It is again for the HG00438 diploid assembly:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_fit_model.png" width="700" height="275">

### 4. Extracting Blocks 

Now we have to assign each coverage value to one of the four components (error, duplicated, haploid, and collapsed). To do so, for each coverage value,
we pick the component with the highest probability. For example for the coverage value, 0, the error component is being picked most of the times (the red line).
In the figure below the coverage intervals are colored based on their assigned component.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_fit_model_colored.png" width="700" height="275">

After assigning the coverage values we then assign the bases of the corresponding assembly to the most probable component. Finally we will have 4 bed 
files each of which points to the regions assinged to a single component.
By passing a coverage file (suffix `.cov`) and its corresponding table file (suffix `.table`) to the binary executable `find_blocks` we can produce those
4 bed files.
````
./find_blocks -c read_alignment.cov -t read_alignment.table -p ${prefix}
````

It outputs four bed files:
````
${prefix}.error.bed
${prefix}.duplicated.bed
${prefix}.haploid.bed
${prefix}.collapsed.bed
````
### 5. Contig-Specific Coverage Analysis

In step 3 we fit a single model for the whole diploid assembly. In step 4 we used that model to partition the assembly into 4 main components. It has been noticed that the model components may change for different regions and it may affect the accuracy of the partitioning process. In order to make the coverage thresholds more sensitive to the local patterns we fit a separate model for each contig. First we split the whole-genome coverage file produced in step 3 into multiple coverage files one for each contig.
```
./split_contigs_cov -c read_alignment.cov -p ${PREFIX}
```
It will produce a list of coverage files like below.
```
HG00438.diploid.hifi.HG00438#1#JAHBCB010000001.1.cov  
HG00438.diploid.hifi.HG00438#1#JAHBCB010000002.1.cov
HG00438.diploid.hifi.HG00438#1#JAHBCB010000003.1.cov
...
HG00438.diploid.hifi.HG00438#2#JAHBCA010000257.1.cov
HG00438.diploid.hifi.HG00438#2#JAHBCA010000258.1#MT.cov
```
To show how the contig-specific models may be different from the whole-genome model, here is an example of a false duplication in the maternal assembly of HG01175 that couldn't be detected using the whole-genome model.
The coverage distributions of three contigs along with the whole assembly are shown below. Two of those contigs are maternal and the other one is paternal. They are containing equivalent regions in the genome but one of the maternal contigs (`2#JAHBCE010000091.1` the red one ) is having a few mega bases of a false duplication of the paternal contig (`1#JAHBCF010000027.1` the yellow one).

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_dip_hifi_contig_cov_dist.png" width="700" height="275">

After fitting the model the duplicated components can reveal such false duplications as you can see in the figure below. `2#JAHBCE010000105.1` does not have any duplicated component but the two others do. Since we know the maternal haplotype has a correct copy of this region in `2#JAHBCE010000105.1` we can deduce that `1#JAHBCF010000027.1` is assembled correctly and `2#JAHBCE010000091.1` is actually containing the false duplication. Without having more information it is not always easy to conclude which copy is the real one and which copy is the false one especially for segmental duplications within a single haplotype.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_contig_fit_model.png" width="700" height="400">

### 6. Combining The Bed Files
One important observation is that for short contigs we don't have a smooth coverage distribution and it is not possible to fit the mixture model. To address this issue we have done the contig-based coverage analysis only for the contigs longer than 5Mb and for the shorter contigs we use the results of the whole-genome-based analysis described previously. So the final results for the current release is a combination of contig- and whole-genome-based coverage anslysis. (Look at the results section, the combined bed files are available in the `combined` subdirectory)


### 7. Correcting The Bed Files Pointing to The Duplicated Regions
1. In some cases the duplicated component is mixed up with the haploid one. It usually happens when the coverage in the haploid component drops systematically and the contig has long stretches of false duplication. One extreme case happened in `HG002#1#JAHKSE010000028.1`. Nearly half of this contig is falsely duplicated but the current results report almost the whole contig as `duplicated`. (https://s3-us-west-2.amazonaws.com/human-pangenomics/submissions/e9ad8022-1b30-11ec-ab04-0a13c5208311--COVERAGE_ANALYSIS_Y1_GENBANK/V1/HiFi/HG002/HG002.diploid.f1_assembly_v2_genbank.hifi.winnowmap_v2.03.cov_dist.pdf). One other indicator of a false duplication is the accumulation of alignments with low MAPQ. We have produced another coverage file only for MAPQ > 20 alignments and intersect it with the current results to correct them. Whenever we see a region flagged as duplicated but it also has more than 5 high quality alignments we flag that region as haploid instead. (Look at the results section, the dup corrected bed files are available in the `dup_corrected` subdirectory)

````
#Get blocks with more than 5 high quality alignments
cat high_mapq.cov | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40)} else if($3 >= 5) {print contig"\t"$1-1"\t"$2}}' | \
            bedtools merge -i - > high_mapq.bed

# do the correction
bedtools subtract -a ${PREFIX}.combined.duplicated.bed -b high_mapq.bed > ${PREFIX}.dup_corrected.duplicated.bed
bedtools intersect -a ${PREFIX}.combined.duplicated.bed -b high_mapq.bed > dup_to_hap.bed
cat dup_to_hap.bed ${PREFIX}.combined.haploid.bed | bedtools sort -i - | bedtools merge -i - > ${PREFIX}.dup_corrected.haploid.bed
````

In the final bed files we have a noticeable number of very short blocks (a few bases or a few tens of bases long). They are very prone to be falsely categorized into one of the components. To increase the specificity we have merged the blocks closer than 100 and then removed the ones shorter than 1Kb.  (Look at the results section, the filtered bed files are available in the `filtered` subdirectory). These numbers may change in the next releases to make the blocks more contiguous.

### Known issues

1. Alignments with low mapping quality are usually happening in the regions with low heterozygosity. The reads with such alignments have to be phased more accurately. Removing (or at least modeling) the read errors and detecting the marker snps are the main steps for finding the correct haplotype of each read, which will be explored for the next releases.
 
2. Some regions are falsely flagged as collapsed. The reason is that the equivalent region in the other haplotype is not assembled correctly so the reads from two haplotypes are aligned to only one of them. This flagging can be useful since it points to a region whose counterpart in the other haplotype is not assembled correctly or not assembled at all. 

    One approach is to correct the coverage in the correctly assembled region. The coverage can be corrected by detecting the marker snps and removing the reads from the wrong halpotype or segment. This is going to be incorporated in the next releases of this analysis. Here is an example of a region with ~40X coverage but after detecting the marker snps (by variant calling) and removing the wrong alignments the coverage has decreased to ~17X which is much closer to the expected coverage (~20X).

   <img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/coverage_correction.png" width="700" height="275">


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

### Acknowledgements

Thanks to Jordan M Eizenga for sharing the source code of EM that fits the mixture model and his useful support for doing this analysis.

