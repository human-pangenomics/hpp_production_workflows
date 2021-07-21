## A Mixture Model-Based Coverage Analysis For HPRC Y1 Assemblies

### Overview
In this page we are explaining the steps toward the coverage analysis of HPRC-Y1 assemblies. Here is a summary of this analysis 
- Align the long reads to each assembly
- Calculate the read coverage of each base of the assembly 
- Fit a mixture model to the coverage distribution
- Extract the blocks that are assigned to the 4 main components of the model; erroreous, haploid, diploid and collapsed.

To be more specific each assembly is partitioned into two (for female samples) or four (for male samples) sets of regions and the third and fourth 
steps above are applied on each set separately. We expect each partition to have different estimated paramters after fitting the model.
There are two main factors that led to this partitioning. The first reason is that the centromeric regions are highly enriched 
in false duplication and collapsing, which is not the case for non-centromeric regions. The other factor is that for male samples
we expect to have half of the sequencing coverage in X and Y chromosomes and using the same model that was fit to the autosomes leads
to overestimating the low-coverage blocks in sex chromosomes.

Each haploid assembly of a female sample is partitioned into 2 sets of regions:
1. non-centromeric (sex + autosome)
2. centromeric (sex + autosome)

Each haploid assembly of a male sample is partitioned into 4 sets of regions:

1. Non-centromeric autosomes + pseudo-atousomal region (PAR)
2. Centromeric autosomes
3. Non-centromeric sex 
4. Centromeric sex

### 1. Read Alignment
The ONT and HiFi reads are aligned to each assembly with winnowmap v2.0.
Here are the main commands used for producing the alignments (adopted from the [winnowmap docs](https://github.com/marbl/Winnowmap)):
```` 
  # making the k-mer table with meryl
  meryl count k=15 output merylDB asm.fa
  meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt
  
  # alignment with winnowmap
  winnowmap -W repetitive_k15.txt -ax [map-ont | map-pb] asm.fa reads.fq.gz | \
    samtools view -hb > read_alignment.bam
````


### 2. Calculating Depth of Coverage
With `samtools depth` we could calculate the the depth of coverage for each base of an assembly. The output of `samtools depth -aa` is like below. (`-aa` 
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
where the name of each contig appears only once before the first block of that contig and the consecutive bases with the same coverage take only one line.
The number that comes after the name of the contig is the contig size. We added the suffix `.cov` to the files with this format.
Here are the main commands for producing the coverage files:
````
samtools depth -aa -Q 20 alignment.bam > alignment.depth
./depth2cov -d alignment.depth -f asm.fa.fai -o read_alignment.cov
````
Note that here the read alignments with mapping quality lower than 20 are filtered out. `-Q 20`
`depth2cov` is a binary executable that converts the output of `samtools depth` to `.cov` format. Its source code is written in C.

### 3. Assembly Alignment and Partitioning


In order to partition an assembly into sex and non-sex or centromeric and non-centromeric regions, 
we align the contigs of that assembly to the chm13v1.0 (+ hg38-Y) reference by winnowmap v2.0. Here are the main commands:
````
# making the k-mer table with meryl
meryl count k=19 output merylDB chm13v1.0_plus_hg38Y.fa
meryl print greater-than distinct=0.9998 merylDB > repetitive_k19.txt

# alignment with winnowmap
winnowmap -W repetitive_k19.txt -ax asm5 chm13v1.0_plus_hg38Y.fa asm.fa | \
  samtools view -hb > asm_alignment.bam
````
Note that the k-mer size for the contig alignment is 19 whereas it was 15 for long read alignment.
After producing the bam file we convert it to the PAF format by `paftools.js` which is available [here](https://github.com/lh3/minimap2/tree/master/misc). Then we pass the paf file
to a python script called `project_blocks.py`. It uses the cigar format in the paf file to project a set of blocks in the reference 
onto the assembly; for this aim we should use `--mode ref2asm`. The desired reference blocks are given as a bed file. Here the bed file may point to 
a specific region in the reference e.g. sex chromosomes. For male samples we prepared 4 distint bed files and for female samples we prepared 2. 
The regions where those bed files are pointing to, have already been mentioned in the overview.
Here is an example:
````
python3 project_blocks.py --mode 'ref2asm' --paf asm_alignment.paf --blocks ref_blocks.bed --outputProjectable projectable.bed --outputProjection asm_blocks.bed
````
For example `ref_blocks.bed` may point to the non-centromeric segments of the sex chromosomes so the `asm_blocks.bed` will point to the assembly regions
mapped to those. That will be used for fitting a separate model to the coverage distribution of these regions.  

### 4. Partitioning The Coverage Files


After the 2nd section, for each assembly and available platform (ONT or HiFi) we will have a coverage file with the suffix `.cov`.
After the 3rd section, for each assembly we will have a partitioning that covers almost all of the contigs.( The contigs with no alignment 
to the reference are excluded from the analysis in the current implementation). A partitioning is represented in the form of 
a set of bed files with no mutual overlap. Each coverage file will be partitioned into 4 (for male) or 2 (for female) smaller coverage files.
For this aim we use a combination of `bedtools` and `awk`:
````
cat read_alignment.cov | \
            awk '{if(substr($1,1,1) == ">") {contig=substr($1,2,40); len_contig=$2} else {print contig"\t"$1-1"\t"$2"\t"$3"\t"len_contig}}' | \
            bedtools intersect -a - -b asm_blocks.bed | \
            awk '{if(contig != $1){contig=$1; print ">"contig"\t"$5}; print $2+1"\t"$3"\t"$4}' > subset_read_alignment.cov
````
For each male assembly we call this command 4 times; once for each bed file. Similarly we call it 2 times for a female assembly.
The suffix of a partitioned file indicates the regions it includes. 
As an example here is a list of the coverage files for two assemblies after being partitioned.
````
HG01258 (A male sample):

HG01258.maternal.ont.autosome_cntr.cov
HG01258.maternal.ont.autosome_nonCntr.cov
HG01258.maternal.ont.sex_cntr.cov
HG01258.maternal.ont.sex_nonCntr.cov
````

````
HG03540 (A female sample):

HG03540.paternal.hifi.diploid_cntr.cov
HG03540.paternal.hifi.diploid_nonCntr.cov
````
### 5. Coverage Distribution and Fitting The Model


After finishing the fourth step we should have 4 and 2 coverage files for each male and female assembly respectively. The frequencies of coverages can
be calculated w/ `cov2counts` which is a binary executable. The source code is written in C.
````
./cov2counts -i read_alignment.cov -o read_alignment.counts
````
The output file with `.counts` suffix is a 2-cloumn tab-delimited file; the first column shows a coverage and the second column shows the frequency of
that coverage in the corresponding partition of the assembly. Using `.counts` files we can produce distribution plots easily. For example below we are 
showing the coverage distributions for the 4 partitions of HG00438 paternal assembly.
<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438.paternal.HiFi.dist.png" width="700" height="375">

Here is a closer look on the lower distributions:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438.paternal.HiFi.dist.zoomed.png" width="700" height="375">

The differences between the 4 distributions clearly show why we did the partitioning in the previous step.
Then we pass each `.counts` file to `fit_model_extra.py` whose job is to fit the mixture model and find the best parameters through 
[Expectation Maximization (EM)](https://en.wikipedia.org/wiki/Expectation%E2%80%93maximization_algorithm). This mixture model
consists of 4 main components and each component represents a specific type of regions.

The 4 components:

1. **Error component** which is modeled by a truncated exponential distribtuion. It represents the regions with very low read support.
2. **Haploid component** which is modeled by a gaussian distribution whose mean is constrained to be half of the diploid component's mean. It should mainly represents the haploid
regions. We expect the false duplicated blocks to also appear in this component. The weight of this component is usually less than 0.01 in non-centromeric partitions but it becomes noticeable
in centromeric ones. It is worth noting that according to the recent 
[T2T paper, The complete sequence of a human genome,](https://www.biorxiv.org/content/10.1101/2021.05.26.445798v1.abstract) there exist
 some satellite arrays (especially HSAT1) where the ONT and HiFi coverage drops systematically due to bias in sample preparation and sequencing.
 As a result this mode should contain a mix of haploid, duplicated and coverage-biased blocks, which are not easy to be distiguished through the current analysis.
3. **Diploid component** which is modeled by a gaussian distribution. It represents blocks with the coverages that we expect for the homozygous blocks of an error-free assembly.
4. **Collpased component** which is actually a set of components each of which follows a gaussian distribution and their means are constrained to be
multiples of the haploid component's mean.

Here is the command that fits the model:
````
python3 fit_model_extra.py --counts read_alignment.counts --output read_alignment
````

Its output is a file with `.table` suffix. It contains a TAB-delimited table with the 7 fields described below:
|Col|Type  |Description                               |
|:--|:----:|:-----------------------------------------|
|coverage  |int|A depth of coverage                       |
|freq  |float   |The frequency of the coverage value in the first column                     |
|fit  |float   |The frequency value fit to the model    |
|error  |float   |The weight of the error component       |
|haploid  |float  |The weight of the haploid component               |
|diploid  |float|The weight of the diploid component                      |
|collapsed  |float   |The weight of the collapsed component                    |

Here is an example of such a table:
````
#coverage	freq	fit	error	haploid	diploid	collapsed
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

In the figure below we are showing the model components after being fit to the diploid centromeres. It is again for the HG00438 paternal assembly:

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438.paternal.HiFi.diploid_cntr.png" width="700" height="375">

### 6. Extracting Blocks

Now we have to assign each coverage value to one of the four components (error, haploid, diploid and collapsed). To do so for each coverage value
we pick the component with the highest weight. For example for the coverage value, 0, the error component is being picked most of the times (the red line).
In the figure below the coverage intervals are colored based on their assigned component.

<img src="https://github.com/human-pangenomics/hpp_production_workflows/blob/asset/coverage/images/HG00438_fit_colored.png" width="700" height="375">

After assigning the coverage values we then assign the bases of the corresponding assembly to the most probable component. Finally we will have 4 bed 
files each of which points to the regions assinged to a single component.
By passing a coverage file (suffix `.cov`) and its corresponding table file (suffix `.table`) to the binary executable `find_blocks` we can produce the 
4 bed files mentioned above.
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

### Known issues


1. The haploid component may contain the regions with a deletion in only one haplotype (or we can equivalently say an insertion in the other haplotype). It may also contain the falsely duplicated regions. In order to separate these two sets of regions one solution is to do the read alignment to both haplotypes at the same time (Thanks to Heng Li for this idea). In that case we expect the main component to be haploid and if there is a homozygous block that exist in both haplotypes the read coverage should be randomly distributed between those two blocks. In another scenario if that homozygous block is falsely duplicated in both haplotypes the coverage analysis should be able to detect that. (This option is currently being explored.)

2. The more coverage we have the more accurate we will be in estimating the parameters of the model. So the samples that have low coverage may not provide a well-fitted coverage distribution. Here is a list of the samples with (<20X) ONT data:

````
HG01123
HG03579
HG03453
HG02622
HG02572
````

3. The contigs that couldn't be mapped to the reference chm13v1.0+hg38Y are not partitioned so they are not included in this analysis. In future developements we have to find which partition those contigs are coming from.

### Data, Source Code and Workflows Availability

The haplotype-resolved assemblies of the HPRC-Y1 samples and their corresponding data sets are available in

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=working/HPRC/

For more details read this github page:

https://github.com/human-pangenomics/HPP_Year1_Assemblies

Please note that for the current results we have used the v2 (and v2.1) of the assemblies. The results will be updated for the genbank assemblies.

The Python scripts, C source codes and the binary files are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/docker/coverage/scripts

The wdl files that have been used for this analysis are available in

https://github.com/human-pangenomics/hpp_production_workflows/tree/asset/coverage/wdl/tasks

### Results Availability

The results are available in 

https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/91b5a5f6-ded2-11eb-886e-0a13c5208311--COVERAGE_ANALYSIS_MODEL_BASED/

The directory structure of the results for a female sample is as below: (For male samples we have the same structure except that the set of partitions are different)
````
.
├── mat
│   ├── HiFi
│   │   ├── beds
│   │   ├── figures
│   │   └── tables
│   └── Ont
│       ├── beds
│       ├── figures
│       └── tables
└── pat
    ├── HiFi
    │   ├── beds
    │   ├── figures
    │   └── tables
    └── Ont
        ├── beds
        ├── figures
        └── tables
````
As it is obvious for each combination of haplotype (mat or pat) and platform (Ont or HiFi) we have the structure below:
````
.
├── beds
├── figures
└── tables
````
It contains 3 directories:

1. `beds` contains two subdirectories:

````
.
├── gt_1k
│   ├── HG00438.maternal.ont.diploid_cntr.collapsed.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_cntr.duplicated.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_cntr.error.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_cntr.haploid.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_nonCntr.collapsed.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_nonCntr.duplicated.gt_1k.bed
│   ├── HG00438.maternal.ont.diploid_nonCntr.error.gt_1k.bed
│   └── HG00438.maternal.ont.diploid_nonCntr.haploid.gt_1k.bed
└── initial
    ├── HG00438.maternal.ont.diploid_cntr.collapsed.bed
    ├── HG00438.maternal.ont.diploid_cntr.duplicated.bed
    ├── HG00438.maternal.ont.diploid_cntr.error.bed
    ├── HG00438.maternal.ont.diploid_cntr.haploid.bed
    ├── HG00438.maternal.ont.diploid_nonCntr.collapsed.bed
    ├── HG00438.maternal.ont.diploid_nonCntr.duplicated.bed
    ├── HG00438.maternal.ont.diploid_nonCntr.error.bed
    └── HG00438.maternal.ont.diploid_nonCntr.haploid.bed
````

In `initial` each bed file points to the blocks of a specific partition that belongs to a specifc component. So for each combination of partition (`diploid_cntr`
 and `diploid_nonCntr`) and components (`error`, `duplicated`, `haploid`, `collapsed`) we will have a separate bed file. The other subdirectory `gt_1k` contains
 the same bed files but the blocks shorter than 1k are filtered out.

2. `figures` contains two subdirectories:

````
.
├── distribution
│   ├── HG00438.maternal.ont
│   └── HG00438.maternal.ont.zoomed.png
└── model
    ├── HG00438.maternal.ont.diploid_cntr.png
    └── HG00438.maternal.ont.diploid_nonCntr.png
````
`distribution` contains two plots of coverage distributions stratified and colored by partition. The first image gives a broad view of the 
coverage distributions and the one with the suffix `.zoomed` is the same plot but zoomed on the lower distributions.
`model` contains one plot per partition. Each plot shows the coverage distribution of a specifc partition and the fitted distribution stratified 
and colored by the component.

3. `tables` contains one `.table` file per partition:
````
.
├── HG00438.maternal.ont.diploid_cntr.table
└── HG00438.maternal.ont.diploid_nonCntr.table
````
We have already explained what each table file contains in section 5.

By going to the github directory `index` you can find the index files that contain the URLs of the results:
````
.
├── female
│   ├── beds
│   │   ├── female.diploid_cntr.gt_1k_bed.index
│   │   ├── female.diploid_cntr.initial_bed.index
│   │   ├── female.diploid_nonCntr.gt_1k_bed.index
│   │   └── female.diploid_nonCntr.initial_bed.index
│   ├── figures
│   │   └── female.figures.index
│   └── tables
│       └── female.tables.index
└── male
    ├── beds
    │   ├── male.autosome_cntr.gt_1k_bed.index
    │   ├── male.autosome_cntr.initial_bed.index
    │   ├── male.autosome_nonCntr.gt_1k_bed.index
    │   ├── male.autosome_nonCntr.initial_bed.index
    │   ├── male.sex_cntr.gt_1k_bed.index
    │   ├── male.sex_cntr.initial_bed.index
    │   ├── male.sex_nonCntr.gt_1k_bed.index
    │   └── male.sex_nonCntr.initial_bed.index
    ├── figures
    │   └── male.figures.index
    └── tables
        └── male.tables.index
````
### Acknowledgements

Thanks to Jordan M Eizenga for sharing the source code of EM that fits the mixture model and his useful support for doing this analysis.
